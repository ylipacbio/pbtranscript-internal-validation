#!/usr/bin/env python

"""Access data files in SMRTLink IsoSeq job"""

import logging
import os.path as op
import json
import re
from collections import defaultdict
from pbcore.io import ContigSet, FastaWriter, FastqWriter
from pbcore.util.Process import backticks
from pbtranscript.Utils import execute
from pbtranscript.io import BLASRM4Reader

log = logging.getLogger(__file__)

__all__ = [
    "consolidate_xml",
    "json_to_attr_dict",
    "get_job_id_from_url",
    "get_host_from_url",
    "get_job_path_from_url",
    "get_subread_xml_from_job_path",
    "reseq",
    "m42coverage",
    "subset_dict",
    "coverage2str"
]


def consolidate_xml(src, dst):
    """Convert input dataset to output fasta|fastq"""
    w = None
    if dst.endswith(".fa") or dst.endswith(".fasta"):
        w = FastaWriter(dst)
        for r in ContigSet(src):
            w.writeRecord(r)
        w.close()
    elif dst.endswith(".fq") or dst.endswith(".fastq"):
        w = FastqWriter(dst)
        for r in ContigSet(src):
            w.writeRecord(r)
        w.close()
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
    host=get_host_from_url(url_path)
    port=8081
    job_id = get_job_id_from_url(url_path)

    cmd="pbservice get-job --host {host} --port {port} {job_id} | grep 'path' |cut -f 2 -d ':' | tr -d ' '".format(host=host, port=port, job_id=job_id)
    ret, code, errmsg = backticks(cmd)
    if code != 0:
        raise RuntimeError("{cmd} failed: {err}".format(cmd=cmd, err=errmsg))
    job_path = ret[0]
    return job_path

def get_subread_xml_from_job_path(job_path):
    """Given a smrtlink job path, return its asscoiated subread xml, either from workflow/datastore.json or pbscala-job.sh"""
    datastore_json_fn = op.join(job_path, 'workflow', 'datastore.json')
    d = {r['fileTypeId']: r['path'] for r in json.load(open(datastore_json_fn, 'r'))['files']}
    if 'PacBio.DataSet.SubreadSet' in d:
        return str(d['PacBio.DataSet.SubreadSet'])
    else: # fall back to get eid_subread from pbscala-job.sh
        try:
            pbscala_sh_fn = op.join(job_path, 'pbscala-job.sh')
            content = ''.join(x for x in open(pbscala_sh_fn, 'r'))
            key = 'eid_subread:'
            if key in content:
                idx = content.find(key)
                return content[idx+len(key):].split(' ')[0].strip()
        except Exception as e:
            raise ValueError("Could not find subread xml from smrtlink job path %s" % job_path)

def reseq(fa_fn, ref_fn, out_m4, nproc=16):
    """Resequencing, output m4"""
    cmd = 'blasr %s %s --out %s --bestn 1 --minMatch 12 --maxMatch 15 --nCandidates 5 --nproc %s -m 4 --header ' % (fa_fn, ref_fn, out_m4, nproc)
    execute(cmd)

def m42coverage(in_m4):
    """Return coverage of each target in BLASR M4 file"""
    coverage_d = defaultdict(lambda: 0)
    for rec in BLASRM4Reader(in_m4):
        coverage_d[rec.sID] += 1
    return coverage_d

def subset_dict(d, selected_keys):
    """Given a dict {key: int} and a list of keys, return subset dict {key: int} for only key in keys.
    If a selected key is not in d, set its value to 0"""
    return {key:(0 if key not in d else d[key]) for key in selected_keys}

def coverage2str(coverage_d):
    """Convert a coverage dict {transcript: coverage} to string."""
    header = 'transcript\tcoverage'
    lines = ['\t'.join([str(transcript), str(coverage)]) for transcript, coverage in coverage_d.items()]
    return '\n'.join([header] + lines)
