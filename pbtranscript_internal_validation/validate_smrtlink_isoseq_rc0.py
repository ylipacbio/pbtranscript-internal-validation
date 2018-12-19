#!/usr/bin/env python

"""
Validate a smrtlink RC0 isoseq job with reference and gencode.
"""
import os
import os.path as op
import argparse
import sys
import logging

from collections import defaultdict
from csv import DictReader
from pbcore.io import FastaReader

from isocollapse.independent.system import execute, realpath, mkdir, rmpath
from isocollapse.libs import AlignmentFile
from pbtranscript.io import ChainConfig, ContigSetReaderWrapper

from pbtranscript.counting.chain_samples import chain_samples
from .ValidationFiles import ValidationFiles, ValidationRunner
from .Constants import ValidationConstants as C
from .Utils import filter_sam_by_targets, map_to_reference
from .io.SMRTLinkIsoSeq3Files import SMRTLinkIsoSeq3Files

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


def validate_with_Gencode(sorted_rep_bam, gencode_gtf, match_out):
    """
    Input:
      sorted_rep_bam -- sorted pbmm2 BAM output mapping (collapsed) representitve isoforms to reference
      eval_dir -- evaluation directory
    Run matchAnnot to compare sorted_rep_bam with gencode v25 and output to eval_dir
    """
    # Convert bam to sam
    sorted_rep_sam = sorted_rep_bam + '.sam'
    log.info("Converting BAM {} to SAM {}".format(sorted_rep_bam, sorted_rep_sam))
    c0 = 'samtools view -h {} -o {}'.format(sorted_rep_bam, sorted_rep_sam)
    execute(c0)

    # MUST remove sorted_rep_bam alignments which map to chromosomes not in gencode gtf
    # (e.g., human transcripts which map to scaffolds rather than chromosome 1-22 and XYM)
    gtf_chrs = get_gtf_chromosomes(gencode_gtf)
    # Must filter alignments mapping to human scaffold, otherwise, matchAnnot would fail
    log.info("Filtering alignments not mapping to %r", gtf_chrs)
    out_sam = sorted_rep_bam + ".gtf_chrs.sam"
    filtered_sam = sorted_rep_bam + ".not_gtf_chrs.sam"
    filter_sam_by_targets(in_sam=sorted_rep_bam, targets=gtf_chrs,
                          out_sam=out_sam, filtered_sam=filtered_sam)

    log.info("Writing matchAnnot output to %s", match_out)
    cmd = "matchAnnot.py --gtf={0} {1} > {2}".format(gencode_gtf, sorted_rep_bam, match_out)
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
    ns = [0, 0, 0, 0, 0, 0]  # ns[i] number of isoforms scores == i
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
    """Sanity check inputs and outputs
    Input SMRTLink job must be an IsoSeq3 job.
    """
    args.smrtlink_job_dir = realpath(args.smrtlink_job_dir)
    is_isoseq3 = op.exists(op.join(args.smrtlink_job_dir, 'tasks', 'isoseq3.tasks.cluster-0'))
    if not is_isoseq3:
        raise ValueError("Job {!r} is not SMRTLink IsoSeq3 job!".format(args.smrtlink_job_dir))
    else:
        log.info("Successfully validated IsoSeq3 job {!r}.".format(args.smrtlink_job_dir))

    args.val_dir = realpath(args.val_dir)

    if not op.exists(args.smrtlink_job_dir):
        raise IOError("SMRTLink job directory %s does not exist" % args.smrtlink_job_dir)

    if not op.exists(args.human_reference):
        raise IOError("Human (SIRV) reference %s/%s does not exist." % (args.human_reference))

    if not op.exists(args.gencode_gtf):
        raise IOError("Gencode gtf file %s does not exist." % args.gencode_gtf)

    log.info("Making validation output directory %s", args.val_dir)
    mkdir(args.val_dir)
    return args


def run(smrtlink_job_dir, val_dir,
        gencode_gtf, human_reference, hg_transcripts_fa, selected_hg_transcripts,
        sirv_truth_dir, sirv_reference, sirv_transcripts_fa,
        sample_name, nproc, make_readlength, args):
    """
    --- args for logging and tracing back.
    1) Collapse HQ isoforms from SMRTLink Iso-Seq3 (w/wo genome) job to hg38.
    """
    runner = ValidationRunner(val_dir, smrtlink_job_dir)

    # make data files and reports to validation dir
    log.info("Making links of smrtlink isoseq outputs")
    runner.make_all_files_from_SMRTLink_job(make_readlength=make_readlength)

    # for human isoforms
    log.info("Compare with gencode annotations.")
    runner.ln_gencode_gtf(gencode_gtf)
    # Collapse HQ isoforms to human and validate with MatchAnot,
    collapse_to_reference(runner, human_reference, runner.collapse_to_hg_dir, nproc)
    log.info("Reseq to human reference.")
    map_to_reference(runner.collapsed_to_hg_rep_fq, human_reference, runner.collapsed_to_hg_rep_sam, nproc)
    log.info("Reseq to human transcripts.")
    runner.reseq_to_human(target_fa=hg_transcripts_fa, selected_transcripts=selected_hg_transcripts.split(','))
    log.info("make readlength plot for Human isoforms")
    runner.make_readlength_csv_for_hg_isoforms()

    validate_human_using_matchAnnot(
            sorted_rep_bam=runner.collapsed_to_hg_rep_sam, gencode_gtf=gencode_gtf,
            matchAnnot_out=runner.matchAnnot_out, hg_report_txt=runner.matchAnnot_out,
            validation_report_csv=runner.validation_report_csv,
            collapsed_to_hg_rep_fq=runner.collapsed_to_hg_rep_fq)

    # for sirv isoforms
    log.info("Reseq to SIRV transcripts.")
    runner.ln_sirv_truth_dir(sirv_truth_dir)
    validate_sirv_isoforms(runner, sirv_reference, sirv_transcripts_fa, sample_name, nproc)
    runner.make_readlength_csv_for_sirv_isoforms()

    # write args and data files to README
    log.info("Writing args and data files to %s", runner.readme_txt)
    runner.make_readme_txt(args=args)


def validate_human_using_matchAnnot(sorted_rep_bam, gencode_gtf, matchAnnot_out, hg_report_txt, validation_report_csv, collapsed_to_hg_rep_fq):
    # Run matchAnnot.py to compare sorted_rep_bam against gencode gtf.
    log.info("Running matchAnnot.py")
    validate_with_Gencode(sorted_rep_bam=sorted_rep_bam, gencode_gtf=gencode_gtf, match_out=matchAnnot_out)

    # Check from matchAnnot.txt % of isoforms with scores 4 or 5
    log.info("Reading matchAnnot reports")
    total_n, ns = check_matchAnnot_out(match_out=matchAnnot_out)

    # collpased isoforms with HIGH matchAnnot score (4 or 5) ? min_percentage
    # ok = (ns[5] + ns[4] >= total_n * min_percentage / 100.0)
    msg = "%s out of %s collapsed isoforms have HIGH MatchAnnot score (>=4)" % (ns[5]+ns[4], total_n)
    log.info(msg)

    # write human related validation metrics to report csv
    writer = open(hg_report_txt, 'w')
    csv_writer = open(validation_report_csv, 'a')
    writer.write(msg + "\n\nDetails:\n")
    csv_writer.write("%s\t%s\n" % ("collapse_to_human.num_isoforms", _get_num_reads(collapsed_to_hg_rep_fq)))
    csv_writer.write("%s\t%s\n" % ("collapse_to_human.num_isoforms_score_ge_4", ns[4] + ns[5]))
    csv_writer.write("%s\t%s\n" % ("collapse_to_human.percentage_isoforms_score_ge_4", (ns[4] + ns[5])/(1.*total_n)))
    for score in range(0, 6):
        writer.write("%s out of %s collapsed isoforms with MatchAnnot score %s\n" % (ns[score], total_n, score))
        csv_writer.write("%s\t%s\n" % ("collapse_to_human.num_isoforms_score_eq_%s" % score, ns[score]))
    writer.close()
    csv_writer.close()


def collapse_to_reference(runner, referenceset_or_fasta, out_dir, nproc):
    """Collapse HQ isoforms to SIRV, get collapsed isoforms in fq and SAM.
    Compare collapsed isoforms against SIRV ground truth,
    return TP, FP, FN
    """
    mkdir(out_dir)

    slfs = runner.sl_job
    aligned_hq_isoforms_bam = op.join(out_dir, 'collapse_hq_transcripts_to_ref.bam')

    # Collapse HQ isoforms fastq to SIRV and make representive isoforms, then map
    # representative isoforms to reference, and sort output BAM
    log.info("Map HQ isoforms to reference {}.".format(referenceset_or_fasta))
    # Collapse HQ isoforms fastq to SIRV
    map_to_reference(readset_or_bam=slfs.hq_transcript_ds,
                     referenceset_or_fasta=referenceset_or_fasta,
                     out_bam=aligned_hq_isoforms_bam, nproc=nproc)

    output_prefix = op.join(out_dir, 'touse.')

    def g(s):
        return output_prefix + s

    log.info("Collapsing HQ isoforms to SIRV.")
    from isocollapse.tasks.collapse_mapped_isoforms import run_main as isocollapse_runner
    import isocollapse.independent.Constants as C
    isocollapse_runner(
        i_hq_transcript_ds=slfs.hq_transcript_ds,
        i_transcriptaln_ds=aligned_hq_isoforms_bam,
        i_all_transcript_ds=slfs.transcript_ds,
        i_flnc_bam=slfs.flnc_bam,
        output_prefix=output_prefix,
        out_isoforms=g('fastq'),
        out_gff=g('gff'),
        out_group=g('group.txt'),
        out_abundance=g('abundance.txt'),
        out_read_stat=g('readstat.txt'),
        out_report_json=g('report.json'),
        min_aln_coverage=C.MIN_ALN_COVERAGE.val,
        min_aln_identity=C.MIN_ALN_IDENTITY.val,
        max_fuzzy_junction=C.MAX_FUZZY_JUNCTION.val,
        allow_extra_5exon=C.ALLOW_EXTRA_5EXON)


def summarize_reseq(bam_fn):
    """Summarize pbmm2 bam file, return (num_mapped_reads, num_mapped_reference)"""
    records = [r for r in AlignmentFile(bam_fn, 'rb')
               if (not r.is_unmapped and not r.is_secondary and not r.is_supplementary)]
    n_mapped_reads = len(set([r.query_name for r in records]))
    n_mapped_refs = len(set([r.target_name for r in records]))
    return (n_mapped_reads, n_mapped_refs)


def validate_sirv_isoforms(runner, sirv_reference, sirv_transcripts_fa, sample_name, nproc):
    """Collapse HQ isoforms to SIRV, get collapsed isoforms in fq and SAM.
    Compare collapsed isoforms against SIRV ground truth,
    return TP, FP, FN
    """
    slfs = runner.sl_job
    # Collapse HQ isoforms fastq to SIRV and make representive isoforms, then map
    # representative isoforms to gmap reference, and sort output SAM (sorted_rep_sam).
    log.info("Collapsing HQ isoforms to SIRV, and mapping representative collapsed isoforms to SIRV.")
    # Collapse HQ isoforms fastq to SIRV
    collapse_to_reference(runner, sirv_reference, runner.collapse_to_sirv_dir, nproc)

    # make a chain config
    mkdir(runner.chain_sample_dir)
    cfg = ChainConfig(sample_names=[C.SIRV_NAME, sample_name],
                      sample_paths=[C.SIRV_TRUTH_DIR, runner.collapse_to_sirv_dir],
                      group_fn=op.basename(runner.collapsed_to_sirv_group),
                      gff_fn=op.basename(runner.collapsed_to_sirv_gff),
                      abundance_fn=op.basename(runner.collapsed_to_sirv_abundance))
    log.info("Write chain config")
    cfg.write(runner.chain_sample_config)

    # same as "chain_samples.py sample.config count_fl --fuzzy_junction=5"
    cwd = os.getcwd()
    # MUST run in chain_sample_dir, all_samples.* will be written to chain_sample_dir/
    os.chdir(runner.chain_sample_dir)
    chain_samples(cfg=cfg, field_to_use='count_fl', max_fuzzy_junction=5)
    os.chdir(cwd)

    # comapre with sirv, get n_total, n_fns, n_fps
    n_total, n_fn, n_fp = compare_with_sirv(chained_ids_fn=runner.chained_ids_txt,
                                            sirv_name=C.SIRV_NAME, sample_name=sample_name)

    report = [('TP', n_total), ('FN', n_fn), ('FP', n_fp)]
    with open(runner.sirv_report_txt, 'w') as f:
        for k, v in report:
            f.write("%s: %s\n" % (k, v))

    # append more sirv validation metrics to report csv
    desc_val_tuples = [
        ("collapse_to_sirv.num_isoforms", _get_num_reads(runner.collapsed_to_sirv_rep_fq)),
        ("collapse_to_sirv.num_TruePositive", n_total),
        ("collapse_to_sirv.num_FalseNegative", n_fn),
        ("collapse_to_sirv.num_FalsePositive", n_fp)
    ]

    # Reseq FLNC/HQ/LQ isoforms to SIRV transcripts
    d = {'reseq_to_sirv.hq_isoforms': (runner.hq_isoforms_fa, runner.hq_sirv_bam),
         'reseq_to_sirv.lq_isoforms': (runner.lq_isoforms_fa, runner.lq_sirv_bam),
         'reseq_to_sirv.isoseq_flnc': (runner.isoseq_flnc_fa, runner.isoseq_flnc_sirv_bam)}
    for name in d.keys():
        fa_fn, bam_fn = d[name][0], d[name][1]
        if os.stat(fa_fn).st_size != 0:
            reseq(fa_fn, sirv_transcripts_fa, bam_fn, nproc)
            n_mapped_reads, n_mapped_refs = summarize_reseq(bam_fn)
        else:
            n_mapped_reads, n_mapped_refs = 0, 0
        desc_val_tuples.append(('%s_n_mapped_reads' % name, n_mapped_reads))
        desc_val_tuples.append(('%s_n_mapped_refs' % name, n_mapped_refs))

    with open(runner.validation_report_csv, 'a') as f:
        for desc, val in desc_val_tuples:
            f.write("%s\t%s\n" % (desc, val))


def compare_with_sirv(chained_ids_fn, sirv_name, sample_name):
    """Compare SIRV results in all_samples.chained_ids.txt against SIRV ground truth,
    return TP, FN, FP"""
    tally = defaultdict(lambda: [])  # SIRV --> list of test ids that hit it (can be redundant sometimes due to fuzzy)

    assert os.path.exists(chained_ids_fn)

    FPs = []
    FNs = []
    for r in DictReader(open(chained_ids_fn, 'r'), delimiter='\t'):
        if r[sirv_name] == 'NA':  # is false positive!
            FPs.append(r[sample_name])
        elif r[sample_name] == 'NA':  # is false negative
            FNs.append(r[sirv_name])
        else:
            tally[r[sirv_name]].append(r[sample_name])

    return len(tally), len(FNs), len(FPs)


def add_human_arguments(parser):
    """Human related args"""
    def get_read_names(fa):
        return ','.join([r.name.split()[0] for r in FastaReader(fa)])
    _group = parser.add_argument_group("Human arguments")
    _group.add_argument("--human_reference", type=str, default=C.HUMAN_SIRV_REFERENCE, help="Human + SIRV reference")
    _group.add_argument("--hg_transcripts_fa", type=str, default=C.HUMAN_TRANSCRIPTS_FA,
                        help="Human reference transcripts FASTA")
    _group.add_argument("--selected_hg_transcripts", type=str, default=get_read_names(C.HUMAN_12_TRANSCRIPTS_FA),
                        help="comma delimited human transcripts to analyze")
    helpstr = "Gencode gtf file containing known human transcripts. default %r" % C.GENCODE_GTF
    _group.add_argument("--gencode_gtf", type=str, default=C.GENCODE_GTF, help=helpstr)
    return parser


def add_sirv_arguments(parser):
    """SIRV related args"""
    _group = parser.add_argument_group("SIRV arguments")
    _group.add_argument("--sirv_transcripts_fa", type=str, default=C.SIRV_TRANSCRIPTS_FA,
                        help="SIRV reference transcripts FASTA")
    helpstr = "SIRV ground truth dir with touse.group.txt, touse.gff, touse.abundance.txt, default %r" % C.SIRV_TRUTH_DIR
    _group.add_argument("--sirv_truth_dir", type=str, default=C.SIRV_TRUTH_DIR, help=helpstr)
    helpstr = "SIRV ground truth reference, default %r" % C.SIRV_REFERENCE,
    _group.add_argument("--sirv_reference", type=str, default=C.SIRV_REFERENCE, help=helpstr)
    return parser


def get_parser():
    """Get argument parser."""
    helpstr = "Validate a SMRTLink Iso-Seq3 job on RC0 sample by comparing against Human GenCode Annotation and SIRV ground truth."
    parser = argparse.ArgumentParser(helpstr, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    helpstr = "Smrtlink Iso-Seq3 job directory"
    parser.add_argument("smrtlink_job_dir", help=helpstr)
    helpstr = "Validation out directory"
    parser.add_argument("val_dir", help=helpstr)

    add_human_arguments(parser)
    add_sirv_arguments(parser)

    parser.add_argument('--make_readlength', default=False, action='store_true', help="Make read length csv.")
    parser.add_argument("--sample_name", type=str, default='sample_name', help="Sample name")
    parser.add_argument("--nproc", type=int, default=C.NPROC, help="NProc")
    return parser


def main(args=sys.argv[1:]):
    """main"""
    args = get_parser().parse_args(args)
    make_sane(args)
    run(smrtlink_job_dir=args.smrtlink_job_dir, val_dir=args.val_dir,
        gencode_gtf=args.gencode_gtf, human_reference=args.human_reference,
        hg_transcripts_fa=args.hg_transcripts_fa,
        selected_hg_transcripts=args.selected_hg_transcripts,
        sirv_truth_dir=args.sirv_truth_dir, sirv_reference=args.sirv_reference,
        sirv_transcripts_fa=args.sirv_transcripts_fa,
        sample_name=args.sample_name, nproc=args.nproc,
        make_readlength=args.make_readlength, args=args)


if __name__ == "__main__":
    main()
