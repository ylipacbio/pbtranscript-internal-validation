#!/usr/bin/env python

"""
Validate a smrtlink RC0 isoseq job with reference and gencode.
"""
import os.path as op
import argparse
import sys
import logging

from pbcore.io import FastaReader

from pbtranscript.Utils import execute, realpath, mkdir, rmpath
from pbtranscript.io import ChainConfig, ContigSetReaderWrapper, BLASRM4Reader
#from pbtranscript.counting.chain_samples import chain_samples
from . import ValidationFiles, ValidationRunner
from . import ValidationConstants as C
from .Utils import reseq, filter_sam_by_targets
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


def summarize_reseq(m4_fn):
    """Summarize resequencing m4 output, return (num_mapped_reads, num_mapped_reference)"""
    records = [r for r in BLASRM4Reader(m4_fn)]
    n_mapped_reads = len(records)
    n_mapped_refs = len(set([r.sID for r in records]))
    return (n_mapped_reads, n_mapped_refs)


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

    if not op.exists(args.reference):
        raise IOError("Human (SIRV) reference %s/%s does not exist." % (args.human_reference))

    if not op.exists(args.gencode_gtf):
        raise IOError("Gencode gtf file %s does not exist." % args.gencode_gtf)

    log.info("Making validation output directory %s", args.val_dir)
    mkdir(args.val_dir)
    return args


def run(smrtlink_job_dir, val_dir, human_reference, gencode_gtf, make_readlength, nproc, args):
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
    collapse_to_reference(smrtlink_job_dir, human_reference, runner.collapse_to_hg_dir, nproc)

    runner.make_readlength_csv_for_hg_isoforms()

    #log.info("Reseq to human transcripts.")
    # runner.reseq_to_human(target_fa=hg_transcripts_fa, selected_transcripts=selected_hg_transcripts.split(','))

    # for sirv isoforms
    # runner.ln_sirv_truth_dir(sirv_truth_dir)
    #validate_sirv_isoforms(SMRTLinkIsoSeqFilesCls, post_mapping_to_genome_runner_f)
    # runner.make_readlength_csv_for_sirv_isoforms()

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


def map_to_reference(readset_or_bam, referenceset_or_fasta, out_bam, nproc):
    """
    Map input readset or bam to reference set or fasta, and create output bam.
    """
    cmd = "pbmm2 {ref} {reads} {out_bam} --sort --preset ISOSEQ -j {nproc}".format(
        ref=referenceset_or_fasta, reads=readset_or_bam, out_bam=out_bam, nproc=nproc)
    log.debug(cmd)
    execute(cmd)


def collapse_to_reference(smrtlink_job_dir, referenceset_or_fasta, out_dir, nproc):
    """Collapse HQ isoforms to SIRV, get collapsed isoforms in fq and SAM.
    Compare collapsed isoforms against SIRV ground truth,
    return TP, FP, FN
    """
    mkdir(out_dir)

    # vfs = ValidationFiles(out_dir)
    slfs = SMRTLinkIsoSeq3Files(smrtlink_job_dir)
    aligned_hq_isoforms_bam = op.join(out_dir, 'aligned')

    # Collapse HQ isoforms fastq to SIRV and make representive isoforms, then map
    # representative isoforms to reference, and sort output BAM
    log.info("Map HQ isoforms to SIRV.")
    # Collapse HQ isoforms fastq to SIRV
    map_to_reference(readset_or_bam=slfs.hq_transcript_ds,
                     referenceset_or_fasta=referenceset_or_fasta,
                     out_bam=aligned_hq_isoforms_bam, nproc=nproc)

    output_prefix = op.join(out_dir, 'touse')

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


def get_parser():
    """Get argument parser."""
    helpstr = "Validate a SMRTLink Iso-Seq3 job on RC0 sample by comparing against Human GenCode Annotation and SIRV ground truth."
    parser = argparse.ArgumentParser(helpstr, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    helpstr = "Smrtlink Iso-Seq3 job directory"
    parser.add_argument("smrtlink_job_dir", help=helpstr)

    helpstr = "Validation out directory"
    parser.add_argument("val_dir", help=helpstr)
    parser.add_argument('--make_readlength', default=False, action='store_true', help="Make read length csv.")

    helpstr = "Gencode gtf file containing known human transcripts. default %r" % C.GENCODE_GTF
    parser.add_argument("--gencode_gtf", type=str, default=C.GENCODE_GTF, help=helpstr)

    helpstr = "SIRV ground truth dir with touse.group.txt, touse.gff, touse.abundance.txt, default %r" % C.SIRV_TRUTH_DIR
    parser.add_argument("--sirv_truth_dir", type=str, default=C.SIRV_TRUTH_DIR, help=helpstr)

    parser.add_argument("--sample_name", type=str, default='sample_name', help="Sample name")
    return parser


def main(args=sys.argv[1:]):
    """main"""
    args = get_parser().parse_args(args)
    make_sane(args)
    run(smrtlink_job_dir=args.smrtlink_job_dir, val_dir=args.val_dir, human_reference=args.human_reference,
        gencode_gtf=args.gencode_gtf, make_readlength=args.make_readlength, nproc=args.nproc, args=args)


if __name__ == "__main__":
    main()
