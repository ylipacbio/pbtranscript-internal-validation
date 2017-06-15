#!/usr/bin/env python

"""
Validate a smrtlink RC0 isoseq job with reference and gencode.
"""
import os
import os.path as op
import random
import argparse
import sys
import logging
from collections import defaultdict
from csv import DictReader

from pbcommand.pb_io.common import load_pipeline_chunks_from_json
import pbcoretools.chunking.chunk_utils as CH
from pbcoretools.chunking.gather import get_datum_from_chunks_by_chunk_key
from pbcore.io import FastaReader

from pbtranscript.RunnerUtils import sge_job_runner
from pbtranscript.ClusterOptions import SgeOptions
from pbtranscript.collapsing import concatenate_sam, sort_sam
from pbtranscript.Utils import execute, realpath, mkdir, rmpath
from pbtranscript.io import ChainConfig, ContigSetReaderWrapper, SMRTLinkIsoSeqFiles, BLASRM4Reader
from pbtranscript.tasks.post_mapping_to_genome import add_post_mapping_to_genome_arguments, \
                post_mapping_to_genome_runner
from pbtranscript.counting.chain_samples import chain_samples
from . import ValidationFiles, ValidationRunner
from . import ValidationConstants as C
from .Utils import reseq
from pbtranscript.collapsing.CollapsingUtils import map_isoforms_and_sort


__author__ = 'etseng@pacb.com, yli@pacb.com'

FORMATTER = op.basename(__name__) + ':%(levelname)s:'+'%(message)s'
logging.basicConfig(level=logging.INFO, format=FORMATTER)
log = logging.getLogger(__name__)


def _get_num_reads(fn):
    """Return number of reads in a FASTA/FASTQ/ContigSet file."""
    reader = ContigSetReaderWrapper(fn)
    n = 0
    for dummy_r in reader:
        n += 1
    return n


def scatter_fastq(fq_fn, nchunks, out_dir):
    """split fq_fn into multiple chunks"""
    mkdir(out_dir)
    scatter_json = op.join(out_dir, "scatter_fastq.json")
    CH.write_fastq_chunks_to_file(scatter_json, fq_fn, nchunks, out_dir, 'chunk', 'fastq')
    chunks = load_pipeline_chunks_from_json(scatter_json)
    o_fastq_fns = [str(fn) for fn in get_datum_from_chunks_by_chunk_key(chunks, '$chunk.fastq_id')]
    return o_fastq_fns


def sge_gmap_cmds(fq_fns, gmap_db, gmap_name, gmap_nproc, sge_queue):
    """Make gmap cmds and submit cmds to sge, wait for them to finish"""
    # a list of scripts to run on sge
    script_fns = [fq_fn + ".map_isoforms_to_genome.sh" for fq_fn in fq_fns]
    sam_fns = [fq_fn + ".sam" for fq_fn in fq_fns]

    def cmd(fq_fn, sam_fn):
        """return a cmd string which maps fq_fn to gmap_db/gmap_name and outputs to sam_fn"""
        return "python -m pbtranscript.tasks.map_isoforms_to_genome %s %s --gmap_name %s --gmap_db %s --gmap_nproc %s" % (fq_fn, sam_fn, gmap_name, gmap_db, gmap_nproc)
    cmds = [cmd(fq_fn, sam_fn) for fq_fn, sam_fn in zip(fq_fns, sam_fns)]
    sge_job_runner(cmds_list=cmds, script_files=script_fns, num_threads_per_job=gmap_nproc,
                   sge_opts=SgeOptions(random.randint(0, 100000), sge_queue=sge_queue),
                   wait_timeout=10000, run_timeout=100000,
                   rescue=None, rescue_times=3)
    return sam_fns


def gather_sam(i_sam_fns, o_sam_fn):
    """Merge sam files in i_sam_fns, outputs to o_sam_fn"""
    unsorted_sam_fn = o_sam_fn + ".unsorted.sam"
    # concatenate sam
    concatenate_sam(i_sam_fns, unsorted_sam_fn)
    # then sort
    sort_sam(unsorted_sam_fn, o_sam_fn)
    # remove intermediate file
    rmpath(unsorted_sam_fn)


def sge_map_isoform_to_genome(fq_fn, gmap_db, gmap_name, o_sam_fn, nchunks, work_dir, gmap_nproc, sge_queue):
    """scatter the task of mapping fq_fn to gmap_db/gmap_name
    into nchunks, run scattered tasks on sge, and gather alignments to o_sam_fn
    """
    o_fq_fns = scatter_fastq(fq_fn, nchunks, work_dir)
    log.info("Scatter %s to %s.", fq_fn, o_fq_fns)
    o_sam_fns = sge_gmap_cmds(o_fq_fns, gmap_db, gmap_name, gmap_nproc, sge_queue)
    gather_sam(i_sam_fns=o_sam_fns, o_sam_fn=o_sam_fn)

def summarize_reseq(m4_fn):
    """Summarize resequencing m4 output, return (num_mapped_reads, num_mapped_reference)"""
    records = [r for r in BLASRM4Reader(m4_fn)]
    n_mapped_reads = len(records)
    n_mapped_refs = len(set([r.sID for r in records]))
    return (n_mapped_reads, n_mapped_refs)

def collapse_to_reference(hq_fq, gmap_db, gmap_name, hq_lq_prefix_pickle,
                          out_dir, args, out_rep_fq, out_rep_sam,
                          out_gff, out_abundance, out_group):
    """
    hq_fq --- HQ isoforms in fastq format.
    out_dir --- where to output results.
    min_count --- minimum # of supportive FLNC reads to call an isoform HQ

    First map HQ isoforms against reference ($gmap_db_dir/$gmap_db_name),
    next collapse HQ isoforms to representative isoforms based on mapping,
    and finally map representative isoforms to reference.

    Return (collapsed isoforms FASTQ, sorted SAM output)
    """
    # Map HQ isoforms to GMAP reference genome
    log.info("Mapping HQ isoforms %s to reference %s/%s.", hq_fq, gmap_db, gmap_name)
    mkdir(out_dir)
    hq_sam = op.join(out_dir, "%s.sorted.sam" % op.basename(hq_fq))
    # GMAP jobs are very slow, decide to run GMAP jobs on SGE instead
    if not args.use_sge:
        map_isoforms_and_sort(input_filename=hq_fq, sam_filename=hq_sam,
                              gmap_db_dir=gmap_db, gmap_db_name=gmap_name, gmap_nproc=C.GMAP_NPROC)
    else:
        sge_map_isoform_to_genome(fq_fn=hq_fq, gmap_db=gmap_db, gmap_name=gmap_name,
                                  o_sam_fn=hq_sam, nchunks=args.sge_nodes,
                                  work_dir=op.join(out_dir, 'sge_map_hq_to_genome'),
                                  gmap_nproc=C.GMAP_NPROC, sge_queue=args.sge_queue)

    log.info("Collapsing and filtering HQ isoforms to create representative isoforms.")
    # Post mapping to genome analysis, including
    #   * collapse polished HQ isoform clusters into groups
    #   * count abundance of collapsed isoform groups
    #   * filter collapsed isoforms based on abundance info
    #rep_fq = op.join(out_dir, "%s.rep.fastq" % out_prefix)
    post_mapping_to_genome_runner(in_isoforms=hq_fq, in_sam=hq_sam,
                                  in_pickle=hq_lq_prefix_pickle,
                                  out_isoforms=out_rep_fq,
                                  out_gff=out_gff,
                                  out_abundance=out_abundance,
                                  out_group=out_group, out_read_stat=None,
                                  min_aln_coverage=args.min_aln_coverage,
                                  min_aln_identity=args.min_aln_identity,
                                  min_flnc_coverage=args.min_flnc_coverage,
                                  max_fuzzy_junction=args.max_fuzzy_junction,
                                  allow_extra_5exon=args.allow_extra_5exon,
                                  min_count=args.min_count)

    # Map representitive isoforms to reference
    log.info("Mapping representative isoforms %s to reference %s/%s", out_rep_fq, gmap_db, gmap_name)
    if not args.use_sge:
        map_isoforms_and_sort(input_filename=out_rep_fq, sam_filename=out_rep_sam,
                              gmap_db_dir=gmap_db, gmap_db_name=gmap_name, gmap_nproc=C.GMAP_NPROC)
    else:
        sge_map_isoform_to_genome(fq_fn=out_rep_fq, gmap_db=gmap_db, gmap_name=gmap_name,
                                  o_sam_fn=out_rep_sam, nchunks=args.sge_nodes,
                                  work_dir=op.join(out_dir, 'sge_map_rep_to_genome'),
                                  gmap_nproc=C.GMAP_NPROC, sge_queue=args.sge_queue)
    return out_rep_fq, out_rep_sam


def filter_sam_by_targets(in_sam, targets, out_sam, filtered_sam):
    """ Read alignments in in_sam, copy alignments with targe name in targets to out_sam,
    write other alignments to filtered_sam."""
    log.info("Writing alignments mapping to %r in %s", targets, out_sam)
    log.info("Writing alignments not mapping to %r in %s", targets, filtered_sam)
    out_writer = open(out_sam, 'w')
    filtered_writer = open(filtered_sam, 'w')
    for line in open(in_sam, 'r'):
        if line.startswith('#') or line.startswith('@'):
            # copy header to both out_sam and filtered_sam
            out_writer.write(line)
            filtered_writer.write(line)
        else:
            if len(line.strip()) == 0:
                continue
            else:
                target = line.split('\t')[2].strip()
                if target in targets:
                    out_writer.write(line)
                else:
                    filtered_writer.write(line)
    out_writer.close()
    filtered_writer.close()

def get_gtf_chromosomes(gtf_fn):
    """Return a set of chromosomes in gtf"""
    chrs = set()
    with open(gtf_fn, 'r') as reader:
        for line in reader:
            if line.startswith('#') or len(line.strip()) == 0:
                continue
            else:
                chrs.add(line.split('\t')[0])
    return chrs


def validate_with_Gencode(sorted_rep_sam, gencode_gtf, match_out):
    """
    Input:
      sorted_rep_sam -- sorted SAM output mapping (collapsed) representitve isoforms to reference
      eval_dir -- evaluation directory
    Run matchAnnot to compare sorted_rep_sam with gencode v25 and output to eval_dir
    """
    # MUST remove sorted_rep_sam alignments which map to chromosomes not in gencode gtf
    # (e.g., human transcripts which map to scaffolds rather than chromosome 1-22 and XYM)
    gtf_chrs = get_gtf_chromosomes(gencode_gtf)
    # Must filter alignments mapping to human scaffold, otherwise, matchAnnot would fail
    log.info("Filtering alignments not mapping to %r", gtf_chrs)
    out_sam = sorted_rep_sam + ".gtf_chrs.sam"
    filtered_sam = sorted_rep_sam + ".not_gtf_chrs.sam"
    filter_sam_by_targets(in_sam=sorted_rep_sam, targets=gtf_chrs,
                          out_sam=out_sam, filtered_sam=filtered_sam)

    log.info("Writing matchAnnot output to %s", match_out)
    cmd = "matchAnnot.py --gtf={0} {1} > {2}".format(gencode_gtf, sorted_rep_sam, match_out)
    execute(cmd)


def check_matchAnnot_out(match_out):
    """
    Check from matchAnnot.txt # of isoforms with scores are 4 or 5
    return total_n and ns, where
    total_n is total number of isoforms
    for i in (0,1,2,4,5) ns[i] is the number of isoforms with matchAnnot score i
    """
    lines = [l for l in open(match_out, 'r') if l.startswith("summary:")]
    total_n = 0
    ns = [0, 0, 0, 0, 0, 0] # ns[i] number of isoforms scores == i
    for l in lines:
        if "isoforms read" in l:
            total_n = [int(s) for s in l.split() if s.isdigit()][0]

        for score in [5, 4, 3, 2, 1, 0]:
            if "isoforms scored %s" % score in l:
                ns[score] = [int(s.replace(',', '')) for s in l.split()
                             if s.replace(',', '').isdigit()][0]
    # sum(ns): total number of isoforms hit at least one gene
    # total_n: total number of isoforms read
    return total_n, ns


def make_sane(args):
    """Sanity check inputs and outputs"""
    args.smrtlink_job_dir = realpath(args.smrtlink_job_dir)
    args.val_dir = realpath(args.val_dir)

    if args.gmap_db is None:
        args.gmap_db = realpath(C.GMAP_DB)
        log.warning("Reset GMAP DB to %s", args.gmap_db)

    if args.hg_gmap_name is None:
        args.hg_gmap_name = C.HUMAN_GMAP_NAME
        log.warning("Reset HUMAN GMAP NAME to %s", args.hg_gmap_name)
    if args.sirv_gmap_name is None:
        args.sirv_gmap_name = C.SIRV_GMAP_NAME
        log.warning("Reset SIRV GMAP NAME to %s", args.sirv_gmap_name)

    if not op.exists(args.smrtlink_job_dir):
        raise IOError("SMRTLink job directory %s does not exist" % args.smrtlink_job_dir)

    if not op.exists(op.join(args.gmap_db, args.hg_gmap_name)):
        raise IOError("Human GMAP reference %s/%s does not exist." % (args.gmap_db, args.hg_gmap_name))

    if not op.exists(op.join(args.gmap_db, args.sirv_gmap_name)):
        raise IOError("SIRV GMAP reference %s/%s does not exist." % (args.gmap_db, args.sirv_gmap_name))

    if not op.exists(args.gencode_gtf):
        raise IOError("Gencode gtf file %s does not exist." % args.gencode_gtf)

    log.info("Making validation output directory %s", args.val_dir)
    mkdir(args.val_dir)
    return args


def run(args):
    """
    Collapse HQ isoforms from SMRTLink Iso-Seq (w/wo genome) job to hg38.
    """
    args = make_sane(args)
    runner = ValidationRunner(args.val_dir, args.smrtlink_job_dir)

    # make data files and reports to validation dir
    log.info("Making links of smrtlink isoseq outputs")
    runner.make_all_files_from_SMRTLink_job(make_readlength=args.make_readlength)

    # for human isoforms
    if args.collapse_to_human:
        runner.ln_gencode_gtf(args.gencode_gtf)
        validate_human_isoforms(args)
        runner.make_readlength_csv_for_hg_isoforms()

    if args.reseq_to_human:
        log.info("Reseq to human transcripts.")
        runner.reseq_to_human(target_fa=args.hg_transcripts_fa, selected_transcripts=args.selected_hg_transcripts.split(','))

    # for sirv isoforms
    if not args.no_sirv:
        runner.ln_sirv_truth_dir(args.sirv_truth_dir)
        validate_sirv_isoforms(args)
        runner.make_readlength_csv_for_sirv_isoforms()

    # write args and  data files to README
    log.info("Writing args and data files to %s", runner.readme_txt)
    runner.make_readme_txt(args=args, collapse_to_human=args.collapse_to_human, reseq_to_human=args.reseq_to_human, no_sirv=args.no_sirv)


def validate_human_isoforms(args):
    """Collapse HQ isoforms to human and validate with MatchAnot,
    NO EXPLICT CRITERIA to determine validation of human isoforms PASS or FAIL.
    """
    vfs = ValidationFiles(args.val_dir)
    slfs = SMRTLinkIsoSeqFiles(args.smrtlink_job_dir)
    # Collapse HQ isoforms fastq to human and make representive isoforms, then map
    # representative isoforms to gmap reference, and sort output SAM (sorted_rep_sam).
    log.info("Collapsing HQ isoforms, and mapping representative collapsed isoforms to reference.")
    dummy_rep_hq, sorted_rep_sam = collapse_to_reference(
        hq_fq=vfs.hq_isoforms_fq, gmap_db=args.gmap_db, gmap_name=args.hg_gmap_name,
        hq_lq_prefix_pickle=slfs.hq_lq_prefix_pickle, out_dir=vfs.collapse_to_hg_dir,
        out_rep_fq=vfs.collapsed_to_hg_rep_fq, out_rep_sam=vfs.collapsed_to_hg_rep_sam,
        out_gff=vfs.collapsed_to_hg_gff, out_abundance=vfs.collapsed_to_hg_abundance,
        out_group=vfs.collapsed_to_hg_group, args=args)

    try:
        # Run matchAnnot.py to compare sorted_rep_sam against gencode gtf.
        log.info("Running matchAnnot.py")
        validate_with_Gencode(sorted_rep_sam=sorted_rep_sam, gencode_gtf=args.gencode_gtf,
                              match_out=vfs.matchAnnot_out)

        # Check from matchAnnot.txt % of isoforms with scores 4 or 5
        log.info("Reading matchAnnot reports")
        total_n, ns = check_matchAnnot_out(match_out=vfs.matchAnnot_out)

        # collpased isoforms with HIGH matchAnnot score (4 or 5) ? min_percentage
        # ok = (ns[5] + ns[4] >= total_n * min_percentage / 100.0)
        msg = "%s out of %s collapsed isoforms have HIGH MatchAnnot score (>=4)" % (ns[5]+ns[4], total_n)
        log.info(msg)

        # write human related validation metrics to report csv
        writer = open(vfs.hg_report_txt, 'w')
        csv_writer = open(vfs.validation_report_csv, 'a')
        writer.write(msg + "\n\nDetails:\n")
        csv_writer.write("%s\t%s\n" % ("collapse_to_human.num_isoforms", _get_num_reads(vfs.collapsed_to_hg_rep_fq)))
        csv_writer.write("%s\t%s\n" % ("collapse_to_human.num_isoforms_score_ge_4", ns[4] + ns[5]))
        csv_writer.write("%s\t%s\n" % ("collapse_to_human.percentage_isoforms_score_ge_4", (ns[4] + ns[5])/(1.*total_n)))
        for score in range(0, 6):
            writer.write("%s out of %s collapsed isoforms with MatchAnnot score %s\n" % (ns[score], total_n, score))
            csv_writer.write("%s\t%s\n" % ("collapse_to_human.num_isoforms_score_eq_%s" % score, ns[score]))
        writer.close()
        csv_writer.close()
    except Exception as e:
        log.warning("Could not validate with matchanot!")

def validate_sirv_isoforms(args):
    """Collapse HQ isoforms to SIRV, get collapsed isoforms in fq and SAM.
    Compare collapsed isoforms against SIRV ground truth,
    return TP, FP, FN
    """
    vfs = ValidationFiles(args.val_dir)
    slfs = SMRTLinkIsoSeqFiles(args.smrtlink_job_dir)

    # Collapse HQ isoforms fastq to SIRV and make representive isoforms, then map
    # representative isoforms to gmap reference, and sort output SAM (sorted_rep_sam).
    log.info("Collapsing HQ isoforms to SIRV, and mapping representative collapsed isoforms to SIRV.")
    # Collapse HQ isoforms fastq to SIRV
    collapse_to_reference(
        hq_fq=vfs.hq_isoforms_fq, gmap_db=args.gmap_db, gmap_name=args.sirv_gmap_name,
        hq_lq_prefix_pickle=slfs.hq_lq_prefix_pickle, out_dir=vfs.collapse_to_sirv_dir,
        out_rep_fq=vfs.collapsed_to_sirv_rep_fq, out_rep_sam=vfs.collapsed_to_sirv_rep_sam,
        out_gff=vfs.collapsed_to_sirv_gff, out_abundance=vfs.collapsed_to_sirv_abundance,
        out_group=vfs.collapsed_to_sirv_group, args=args)

    # make a chain config
    mkdir(vfs.chain_sample_dir)
    cfg = ChainConfig(sample_names=[C.SIRV_NAME, args.sample_name],
                      sample_paths=[C.SIRV_TRUTH_DIR, vfs.collapse_to_sirv_dir],
                      group_fn=op.basename(vfs.collapsed_to_sirv_group),
                      gff_fn=op.basename(vfs.collapsed_to_sirv_gff),
                      abundance_fn=op.basename(vfs.collapsed_to_sirv_abundance))
    log.info("Write chain config")
    cfg.write(vfs.chain_sample_config)

    # same as "chain_samples.py sample.config count_fl --fuzzy_junction=5"
    cwd = os.getcwd()
    # MUST run in chain_sample_dir, all_samples.* will be written to chain_sample_dir/
    os.chdir(vfs.chain_sample_dir)
    chain_samples(cfg=cfg, field_to_use='count_fl', max_fuzzy_junction=5)
    os.chdir(cwd)

    # comapre with sirv, get n_total, n_fns, n_fps
    n_total, n_fn, n_fp = compare_with_sirv(chained_ids_fn=vfs.chained_ids_txt,
                                            sirv_name=C.SIRV_NAME, sample_name=args.sample_name)

    report = [('TP', n_total), ('FN', n_fn), ('FP', n_fp)]
    with open(vfs.sirv_report_txt, 'w') as f:
        for k, v in report:
            f.write("%s: %s\n" % (k, v))

    # append more sirv validation metrics to report csv
    desc_val_tuples = [
            ("collapse_to_sirv.num_isoforms", _get_num_reads(vfs.collapsed_to_sirv_rep_fq)),
            ("collapse_to_sirv.num_TruePositive", n_total),
            ("collapse_to_sirv.num_FalseNegative", n_fn),
            ("collapse_to_sirv.num_FalsePositive", n_fp)
    ]

    # Reseq FLNC/HQ/LQ isoforms to SIRV transcripts
    d = {'reseq_to_sirv.hq_isoforms': (vfs.hq_isoforms_fa, vfs.hq_sirv_m4),
         'reseq_to_sirv.lq_isoforms': (vfs.lq_isoforms_fa, vfs.lq_sirv_m4),
         'reseq_to_sirv.isoseq_flnc': (vfs.isoseq_flnc_fa, vfs.isoseq_flnc_sirv_m4)}
    for name in d.keys():
        fa_fn, m4_fn = d[name][0], d[name][1]
        reseq(fa_fn, args.sirv_transcripts_fa, m4_fn)
        n_mapped_reads, n_mapped_refs = summarize_reseq(m4_fn)
        desc_val_tuples.append(('%s_n_mapped_reads' % name, n_mapped_reads))
        desc_val_tuples.append(('%s_n_mapped_refs' % name, n_mapped_refs))

    with open(vfs.validation_report_csv, 'a') as f:
        for desc, val in desc_val_tuples:
            f.write("%s\t%s\n" % (desc, val))


def compare_with_sirv(chained_ids_fn, sirv_name, sample_name):
    """Compare SIRV results in all_samples.chained_ids.txt against SIRV ground truth,
    return TP, FN, FP"""
    tally = defaultdict(lambda: [])# SIRV --> list of test ids that hit it (can be redundant sometimes due to fuzzy)

    assert os.path.exists(chained_ids_fn)

    FPs = []
    FNs = []
    for r in DictReader(open(chained_ids_fn, 'r'), delimiter='\t'):
        if r[sirv_name] == 'NA': # is false positive!
            FPs.append(r[sample_name])
        elif r[sample_name] == 'NA': # is false negative
            FNs.append(r[sirv_name])
        else:
            tally[r[sirv_name]].append(r[sample_name])

    return len(tally), len(FNs), len(FPs)

def add_reseq_arguments(parser):
    """Resequencing HQ/LQ isoforms against SIRV and Gencode transcripts"""
    reseq_group = parser.add_argument_group("Resequencing arguments")
    reseq_group.add_argument("--sirv_transcripts_fa", type=str, default=C.SIRV_TRANSCRIPTS_FA, help="SIRV reference transcripts FASTA")
    reseq_group.add_argument("--hg_transcripts_fa", type=str, default=C.HUMAN_TRANSCRIPTS_FA, help="Human reference transcripts FASTA")
    def get_read_names(fa):
        return ','.join([r.name.split()[0] for r in FastaReader(fa)])
    reseq_group.add_argument("--selected_hg_transcripts", type=str, default=get_read_names(C.HUMAN_12_TRANSCRIPTS_FA), help="comma delimited human transcripts to analyze")
    return parser

def add_gmap_arguments(parser):
    """Get gmap parser"""
    gmap_group = parser.add_argument_group("GMAP arguments")
    gmap_group.add_argument("--gmap_db", type=str, default=C.GMAP_DB, help="GMAP DB Path, default: %r" % C.GMAP_DB)
    gmap_group.add_argument("--hg_gmap_name", "--gmap_name", dest="hg_gmap_name", type=str,
                            default=C.HUMAN_GMAP_NAME, help="Human GMAP DB name, default: %r" % C.HUMAN_GMAP_NAME)
    gmap_group.add_argument("--sirv_gmap_name", type=str, default=C.SIRV_GMAP_NAME,
                            help="SIRV GMAP DB name, default=%r" % C.SIRV_GMAP_NAME)
    gmap_group.add_argument("--gmap_nproc", type=int, default=C.GMAP_NPROC,
                            help="Number of CPUs, default: %s" % C.GMAP_NPROC)
    return parser


def add_sge_arguments(parser):
    """Get sge parser"""
    sge_group = parser.add_argument_group("SGE arguments")
    sge_group.add_argument("--use_sge", default=False, action='store_true', help="Submit jobs to SGE")
    sge_group.add_argument("--sge_queue", type=str, default='default',
                           help="SGE queue to submit GMAP jobs, e.g., 'default', 'production'")
    sge_group.add_argument("--sge_nodes", type=int, default=12,
                           help="Number of sge nodes to use, default: %s" % 12)
    return parser


def get_parser():
    """Get argument parser."""
    helpstr = "Validate an SMRTLink Iso-Seq job of RC0 sample by comparing with Human GenCode Annotation and SIRV ground truth"
    parser = argparse.ArgumentParser(helpstr, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    helpstr = "Smrtlink Iso-Seq job directory"
    parser.add_argument("smrtlink_job_dir", help=helpstr)

    helpstr = "Validation out directory"
    parser.add_argument("val_dir", help=helpstr)

    parser.add_argument('--no_sirv', default=False, action='store_true', help="Don't compare against SIRV.")
    parser.add_argument('--make_readlength', default=False, action='store_true', help="Make read length csv.")
    parser.add_argument('--collapse_to_human', action='store_true', help="Collapse HQ isoforms to human genome. NOTE that it will be SLOW")
    parser.add_argument('--reseq_to_human', action='store_true', help="Reseq FLNC and HQ isoforms to human transcripts.")

    parser = add_post_mapping_to_genome_arguments(parser)
    parser = add_gmap_arguments(parser)
    parser = add_reseq_arguments(parser)
    parser = add_sge_arguments(parser)

    helpstr = "Gencode gtf file containing known human transcripts. default %r" % C.GENCODE_GTF
    parser.add_argument("--gencode_gtf", type=str, default=C.GENCODE_GTF, help=helpstr)

    helpstr = "SIRV ground truth dir with touse.group.txt, touse.gff, touse.abundance.txt, default %r" % C.SIRV_TRUTH_DIR
    parser.add_argument("--sirv_truth_dir", type=str, default=C.SIRV_TRUTH_DIR, help=helpstr)

    parser.add_argument("--sample_name", type=str, default='sample_name', help="Sample name")
    return parser

def main(args=sys.argv[1:]):
    """main"""
    run(get_parser().parse_args(args))

if __name__ == "__main__":
    main()
