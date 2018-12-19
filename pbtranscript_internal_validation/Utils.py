#!/usr/bin/env python

"""Access data files in SMRTLink IsoSeq job"""

import logging
import os
import os.path as op
import json
import re
from collections import defaultdict
from pbcore.io import ContigSet, FastaWriter, FastqWriter
from isocollapse.independent.system import execute, realpath, mv, backticks
from isocollapse.libs import AlignmentFile

log = logging.getLogger(__file__)


def _bam2fastx(exe, src, dst):
    # generate flnc.fasta from flnc.bam
    tmp = dst + 'tmp'
    cmd = 'rm -f {dst} && {exe} {src} -o {tmp} && gunzip {tmp}.fasta.gz && mv {tmp}.fasta {dst}'.format(
        dst=dst, src=src, tmp=tmp, exe=exe)
    execute(cmd)


def consolidate_xml(src, dst):
    """Convert input contigset or subreadset/ccs/transcriptset/alignmentset to fasta|fastq"""
    def _contigset2fastx(cls, src, dst):
        with cls(dst) as writer:
            for r in ContigSet(src):
                writer.writeRecord(r)

    def _convertible_to_bam(s):
        suffices = ['bam', 'alignmentset.xml', 'subreadset.xml', 'consensusreadset.xml',
                    'transcriptset.xml']
        return any([s.endswith(suffix) for suffix in suffices])

    if dst.endswith('.fa') or dst.endswith('.fasta'):
        if _convertible_to_bam(src):
            _bam2fastx('bam2fasta', src, dst)
        else:
            _contigset2fastx(FastaWriter, src, dst)
    elif dst.endswith(".fq") or dst.endswith(".fastq"):
        if _convertible_to_bam(src):
            _bam2fastx('bam2fastq', src, dst)
        else:
            _contigset2fastx(FastqWriter, src, dst)
    else:
        raise ValueError("Output file %s must be either fasta or fastq", dst)


def load_json(json_fn):
    """"Read a json to a dict"""
    with open(json_fn) as json_data:
        return json.load(json_data)


def json_to_attr_dict(json_fn):
    """Convert a json report to a dict {id: value}."""
    ret = dict()
    d = load_json(json_fn)
    try:
        for attr in d['attributes']:
            ret[attr['id']] = attr['value']
    except ValueError:
        log.warning("Warning: could not get attributes from json report %s", json_fn)
    return ret


def get_job_id_from_url(url_path):
    """Given a smrtlink job url (e.g., https://smrtlink-alpha.nanofluidics.com:8243/sl/#/analysis/job/13695), return job id 13695"""
    return url_path.split('/')[-1]


def get_host_from_url(url_path):
    """Given a smrtlink job url, return host (e.g., 'http://smrtlink-alpha')"""
    return re.sub(r'https', 'http', ':'.join(url_path.split(':')[0:2]))


def get_job_path_from_url(url_path):
    """Given a smrtlink job url, return file path to the job, e.g., '/pbi/.../013695'"""
    host = get_host_from_url(url_path)
    port = 8081
    job_id = get_job_id_from_url(url_path)

    cmd = "pbservice get-job --host {host} --port {port} {job_id} | grep 'path' |cut -f 2 -d ':' | tr -d ' '".format(
        host=host, port=port, job_id=job_id)
    ret, code, errmsg = backticks(cmd)
    if code != 0:
        raise RuntimeError("{cmd} failed: {err}".format(cmd=cmd, err=errmsg))
    job_path = ret[0]
    return job_path


def get_subread_xml_from_job_path(job_path):
    """Given a smrtlink job path, return its asscoiated subread xml, either from workflow/datastore.json or pbscala-job.sh"""
    datastore_json_fn = op.join(job_path, 'workflow', 'datastore.json')
    d = {}
    try:
        d = {r['fileTypeId']: r['path'] for r in json.load(open(datastore_json_fn, 'r'))['files']}
    except Exception:
        pass
    if 'PacBio.DataSet.SubreadSet' in d:
        return str(d['PacBio.DataSet.SubreadSet'])
    else:  # fall back to get eid_subread from pbscala-job.sh
        try:
            pbscala_sh_fn = op.join(job_path, 'pbscala-job.sh')
            content = ''.join(x for x in open(pbscala_sh_fn, 'r'))
            key = 'eid_subread:'
            if key in content:
                idx = content.find(key)
                return content[idx+len(key):].split(' ')[0].strip()
        except Exception as e:
            try:
                d = op.join(job_path, 'entry-points')
                assert op.exists(d)
                return [op.join(d, f) for f in os.listdir(d)
                        if op.isfile(op.join(d, f)) and f.endswith('subreadset.xml')][0]
            except Exception as e:
                raise ValueError("Could not find subread xml from smrtlink job path %s" % job_path)


def bam2coverage(in_bam):
    coverage_d = defaultdict(lambda: 0)
    with AlignmentFile(in_bam, 'rb') as reader:
        for r in reader:
            if not r.is_unmapped and not r.is_secondary and not r.is_supplementary:
                coverage_d[r.reference_name] += 1
    return coverage_d


def subset_dict(d, selected_keys):
    """Given a dict {key: int} and a list of keys, return subset dict {key: int} for only key in keys.
    If a selected key is not in d, set its value to 0"""
    return {key: (0 if key not in d else d[key]) for key in selected_keys}


def coverage2str(coverage_d):
    """Convert a coverage dict {transcript: coverage} to string."""
    header = 'transcript\tcoverage'
    lines = ['\t'.join([str(transcript), str(coverage)]) for transcript, coverage in coverage_d.items()]
    return '\n'.join([header] + lines)


def filter_sam_by_targets(in_sam, targets, out_sam, filtered_sam):
    """ Read alignments in in_sam, copy alignments with targe name in targets to out_sam,
    write other alignments to filtered_sam."""
    if in_sam.endswith('.bam'):
        cmd = 'samtools view -h {} -o {}.sam'.format(in_sam, in_sam)
        execute(cmd)
        in_sam = in_sam + '.sam'

    def ref_to_chr_ref(ref):
        if 'chr{}'.format(ref) in targets:
            return 'chr{}'.format(ref)
        return ref

    log.info("Writing alignments mapping to %r in %s", targets, out_sam)
    log.info("Writing alignments not mapping to %r in %s", targets, filtered_sam)
    out_writer = open(out_sam, 'w')
    filtered_writer = open(filtered_sam, 'w')
    for line in open(in_sam, 'r'):
        if line.startswith('#') or line.startswith('@'):
            # copy header to both out_sam and filtered_sam
            if line.startswith('@SQ\tSN'):
                ref = line.split('\t')[1][len('SN:'):]
                ref = ref_to_chr_ref(ref)
                newline = '@SQ\tSN:{}\t'.format(ref) + '\t'.join(line.split('\t')[2:])
                line = newline
            out_writer.write(line)
            filtered_writer.write(line)
        else:
            if len(line.strip()) == 0:
                continue
            else:
                ref = line.split('\t')[2].strip()
                ref = ref_to_chr_ref(ref)
                if ref in targets:
                    out_writer.write(line)
                else:
                    filtered_writer.write(line)
    out_writer.close()
    filtered_writer.close()


def map_to_reference(readset_or_bam, referenceset_or_fasta, out_bam, nproc):
    """
    Map input readset or bam to reference set or fasta, and create output bam.
    """
    cmd = "pbmm2 align {reads} {ref} {out_bam} --sort --preset ISOSEQ -j {nproc}".format(
        ref=referenceset_or_fasta, reads=readset_or_bam, out_bam=out_bam, nproc=nproc)
    log.debug(cmd)
    execute(cmd)

def reseq(readset_or_bam, referenceset_or_fasta, out_bam, nproc):
    """
    pbmm2 align query.fasta ref.fasta out.bam has a bug displaying incorrect SQ lines.
    Have to fallback to minimap2
    """
    def to_fasta(referenceset_or_fasta):
        if referenceset_or_fasta.endswith('.fasta') or referenceset_or_fasta.endswith('.fa'):
            return referenceset_or_fasta
        elif referenceset_or_fasta.endswith('.xml'):
            from pbcore.io import ReferenceSet
            return ReferenceSet(referenceset_or_fasta).toExternalFiles()[0]
        else:
            raise ValueError("{} is not a reference!".format(referenceset_or_fasta))

    ref_fasta = to_fasta(referenceset_or_fasta)
    isoseq_args = '-ax splice -uf -C5 '
    cmd = "minimap2 -t {nproc} {args} -a {ref} {reads} > {out_bam}.sam && samtools view -h {out_bam}.sam -o {out_bam}".format(
        ref=ref_fasta, reads=readset_or_bam, out_bam=out_bam, nproc=nproc, args=isoseq_args)
    log.debug(cmd)
    execute(cmd)
