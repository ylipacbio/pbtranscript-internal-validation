#!/usr/bin/env python

"""Define isoseq validation output files."""
from __future__ import print_function

import logging
import os.path as op
from collections import defaultdict
from pbcore.io import ContigSet, FastaReader, FastqReader
from pbtranscript.io import ContigSetReaderWrapper
from isocollapse.independent.system import realpath, rmpath, lnabs, mkdir, execute, touch
from .Utils import (consolidate_xml, json_to_attr_dict, get_subread_xml_from_job_path,
                    coverage2str, subset_dict, bam2coverage, reseq)
from .io.SMRTLinkIsoSeq3Files import SMRTLinkIsoSeq3Files

FORMATTER = op.basename(__file__) + ':%(levelname)s:'+'%(message)s'
logging.basicConfig(level=logging.INFO, format=FORMATTER)
log = logging.getLogger(__name__)


__all__ = ["ValidationFiles", "ValidationRunner"]


def make_readlength_csv(fasta_fn, csv_fn):
    """Make read length csv file for fasta file"""
    log.info("Making read length csv file %s from %s", csv_fn, fasta_fn)
    rmpath(csv_fn)
    with open(csv_fn, 'w') as writer:
        writer.write("'name'\t'readlength'\n")
        cls = FastaReader if (fasta_fn.endswith('.gz') or fasta_fn.endswith('.fasta')) else ContigSetReaderWrapper
        for read in cls(fasta_fn):
            writer.write('%s\t%s\n' % (read.name.split()[0], len(read.sequence)))


def tolist(x):
    """convert x to a list"""
    return x if isinstance(x, list) else [x]


class ValidationFiles(object):

    """Simple class for all isoseq validation output files."""

    def __init__(self, root_dir):
        self.root_dir = realpath(root_dir)

    @property
    def csv_dir(self):
        """dir for keeping all csv output"""
        return op.join(self.root_dir, 'csv')

    @property
    def fasta_dir(self):
        """dir for keeping all isoseq fasta files"""
        return op.join(self.root_dir, 'fasta')

    @property
    def validation_report_csv(self):
        """file path to validation report"""
        return op.join(self.root_dir, "isoseq_rc0_validation_report.csv")

    @property
    def sirv_report_txt(self):
        """Report for SIRV isoforms"""
        return op.join(self.root_dir, 'SIRV_evaluation_summary.txt')

    @property
    def hg_report_txt(self):
        """Report for hg isoforms"""
        return op.join(self.root_dir, 'human_evaluation_summary.txt')

    @property
    def gencode_gtf(self):
        """file path to human gencode_gtf file"""
        return op.join(self.root_dir, "gencode.annotation.gtf")

    @property
    def sirv_truth_dir(self):
        """dir for keeping sirv ground truth files,
        including touse.gff, touse.ground.txt, touse.count.txt"""
        return op.join(self.root_dir, "SIRV")

    @property
    def readme_txt(self):
        """README file."""
        return op.join(self.root_dir, 'README_.txt')

    @property
    def chain_sample_dir(self):
        """chain sample dir"""
        return op.join(self.root_dir, 'chain_sample')

    @property
    def chain_sample_config(self):
        """chain sample config"""
        return op.join(self.chain_sample_dir, 'chain_sample.contig')

    @property
    def chained_ids_txt(self):
        """all_samples.chained_ids.txt"""
        return op.join(self.chain_sample_dir, 'all_samples.chained_ids.txt')

    @property
    def polymerase_readlength_csv(self):
        """Return polymerase read length csv"""
        return op.join(self.csv_dir, "polymerase_readlength.csv")

    @property
    def ccs_readlength_csv(self):
        """file path to CCS read length csv"""
        return op.join(self.csv_dir, "ccs_readlength.csv")

    @property
    def flnc_readlength_csv(self):
        """file path to FLNC read length csv"""
        return op.join(self.csv_dir, "flnc_readlength.csv")

    @property
    def consensus_isoforms_readlength_csv(self):
        """file path to consensus isoforms"""
        return op.join(self.csv_dir, "consensus_isoforms_readlength.csv")

    @property
    def hq_readlength_csv(self):
        """file path to HQ read length"""
        return op.join(self.csv_dir, "hq_readlength.csv")

    @property
    def lq_readlength_csv(self):
        """file path to LQ isoforms read length csv"""
        return op.join(self.csv_dir, "lq_readlength.csv")

    @property
    def collapsed_isoforms_readlength_csv(self):
        """file path to collapsed isoforms"""
        return op.join(self.csv_dir, "collapsed_isoforms_readlength.csv")

    @property
    def ccs_fa(self):
        """file path to concatenated ccs.fasta"""
        return op.join(self.fasta_dir, "ccs.fasta")

    @property
    def isoseq_flnc_fa(self):
        """file path to concatenated isoseq_flnc.fasta"""
        return op.join(self.fasta_dir, "isoseq_flnc.fasta")

    @property
    def consensus_isoforms_fa(self):
        """file path to consensus_isoforms.fasta"""
        return op.join(self.fasta_dir, "consensus_isoforms.fasta")

    @property
    def reseq_to_sirv_dir(self):
        return op.join(self.root_dir, 'reseq_to_sirv')

    @property
    def hq_sirv_bam(self):
        return op.join(self.reseq_to_sirv_dir, "hq_isoforms.sirv.bam")

    @property
    def lq_sirv_bam(self):
        return op.join(self.reseq_to_sirv_dir, "lq_isoforms.sirv.bam")

    @property
    def isoseq_flnc_sirv_bam(self):
        return op.join(self.reseq_to_sirv_dir, "isoseq_flnc.sirv.bam")

    @property
    def hq_isoforms_fa(self):
        """file path to hq isoforms.fasta"""
        return op.join(self.fasta_dir, "hq_isoforms.fasta")

    @property
    def hq_isoforms_fq(self):
        """file path to hq isoforms.fastq"""
        return op.join(self.fasta_dir, "hq_isoforms.fastq")

    @property
    def lq_isoforms_fa(self):
        """file path to lq isoforms.fasta"""
        return op.join(self.fasta_dir, "lq_isoforms.fasta")

    @property
    def lq_isoforms_fq(self):
        """file path to lq isoforms.fastq"""
        return op.join(self.fasta_dir, "lq_isoforms.fastq")

    @property
    def collapse_to_hg_dir(self):
        """dir for saving files collapsing HQ isoforms to human reference"""
        return op.join(self.root_dir, "collapse_to_hg")

    @property
    def collapse_to_sirv_dir(self):
        """dir for saving files collapsing HQ isoforms to SIRV reference"""
        return op.join(self.root_dir, "collapse_to_sirv")

    @property
    def collapsed_to_hg_rep_fq(self):
        """Representative isoforms collapsed to human genome"""
        return op.join(self.collapse_to_hg_dir, "touse.fastq")

    @property
    def collapsed_to_sirv_rep_fq(self):
        """Representative isoforms collapsed to sirv genome"""
        return op.join(self.collapse_to_sirv_dir, "touse.fastq")

    @property
    def collapsed_to_hg_rep_bam(self):
        """sorted alignments mapping isoforms to SIRV"""
        return op.join(self.collapse_to_hg_dir, "touse.rep.bam")

    @property
    def collapsed_to_hg_gff(self):
        """gff file of rep isoforms collapsed to hg"""
        return op.join(self.collapse_to_hg_dir, "touse.gff")

    @property
    def collapsed_to_sirv_gff(self):
        """gff file of rep isoforms collapsed to sirv"""
        return op.join(self.collapse_to_sirv_dir, "touse.gff")

    @property
    def collapsed_to_hg_abundance(self):
        """abundance file of rep isoforms collapsed to hg"""
        return op.join(self.collapse_to_hg_dir, "touse.abundance.txt")

    @property
    def collapsed_to_sirv_abundance(self):
        """abundance file of rep isoforms collapsed to sirv"""
        return op.join(self.collapse_to_sirv_dir, "touse.abundance.txt")

    @property
    def collapsed_to_hg_group(self):
        """group file of rep isoforms collapsed to hg"""
        return op.join(self.collapse_to_hg_dir, "touse.group.txt")

    @property
    def collapsed_to_sirv_group(self):
        """abundance file of rep isoforms collapsed to sirv"""
        return op.join(self.collapse_to_sirv_dir, "touse.group.txt")

    @property
    def matchAnnot_out(self):
        """MatchAnot output"""
        return self.collapsed_to_hg_rep_bam + ".matchAnnot.txt"

    @property
    def reseq_to_hg_dir(self):
        """Directory saving files resequencing to human transcripts"""
        return op.join(self.root_dir, "reseq_to_hg")

    @property
    def hq_reseq_to_hg_bam(self):
        """Resequence HQ isoforms to human transcripts using blasr, outputs in bam"""
        return op.join(self.reseq_to_hg_dir, "hq_to_hg.bam")

    @property
    def hq_reseq_to_hg_selected_transcripts_csv(self):
        return op.join(self.reseq_to_hg_dir, 'hq_reseq_to_hg_selected_transcripts.csv')

    @property
    def flnc_reseq_to_hg_bam(self):
        """Reseq flnc reads to human transcripts using blasr, outputs in bam"""
        return op.join(self.reseq_to_hg_dir, "flnc_to_hg.bam")

    @property
    def flnc_reseq_to_hg_selected_transcripts_csv(self):
        return op.join(self.reseq_to_hg_dir, 'flnc_reseq_to_hg_selected_transcripts.csv')

    @property
    def collapsed_to_sirv_rep_bam(self):
        """sorted alignments mapping isoforms to SIRV"""
        return op.join(self.collapse_to_sirv_dir, "touse.rep.bam")

    @property
    def collapsed_to_hg_rep_readlength_csv(self):
        """readlength csv of representative isoforms collapsed to hg"""
        return op.join(self.csv_dir, "collapsed_to_hg_rep_readlength.csv")

    @property
    def collapsed_to_sirv_rep_readlength_csv(self):
        """readlength csv of representative isoforms collapsed to hg"""
        return op.join(self.csv_dir, "collapsed_to_sirv_rep_readlength.csv")

    def __str__(self):
        return "IsoSeq validation files under: %s" % self.root_dir


class ValidationRunner(ValidationFiles):
    """Validating from a smrtlink job dir."""

    def __init__(self, root_dir, smrtlink_job_dir):
        super(ValidationRunner, self).__init__(root_dir=root_dir)
        self.smrtlink_job_dir = op.join(smrtlink_job_dir)
        self.sl_job = SMRTLinkIsoSeq3Files(self.smrtlink_job_dir)

    @property
    def all_files(self):
        """Return a list of all files and dirs as [(file_description, file_path)]"""
        return self.common_files + self.collapse_human_files + self.reseq_human_files + self.sirv_files

    @property
    def common_files(self):
        """all common files, not validation files."""
        return [
            ("root_dir", self.root_dir),
            ("smrtlink_dir", self.smrtlink_job_dir),
            ("validation_report_csv", self.validation_report_csv),
            ("polymerase_readlength_csv", self.polymerase_readlength_csv),
            ("ccs_readlength_csv", self.ccs_readlength_csv),
            ("flnc_readlength_csv", self.flnc_readlength_csv),
            ("consensus_isoforms_readlength_csv", self.consensus_isoforms_readlength_csv),
            ("hq_readlength_csv", self.hq_readlength_csv),
            ("lq_readlength_csv", self.lq_readlength_csv),
            ("isoseq_flnc_fasta", self.isoseq_flnc_fa),
            ("consensus_isoforms_fasta", self.consensus_isoforms_fa),
            ("hq_isoforms_fasta", self.hq_isoforms_fa)
        ]

    @property
    def collapse_human_files(self):
        """human related files."""
        return [
            ("gencode_gtf", self.gencode_gtf),
            ("collapsed_to_hg_rep_fastq", self.collapsed_to_hg_rep_fq),
            ("collapsed_to_hg_rep_bam", self.collapsed_to_hg_rep_bam),
            ("collapsed_to_hg_rep_readlength_csv", self.collapsed_to_hg_rep_readlength_csv),
            ("matchAnnot_out", self.matchAnnot_out)
        ]

    @property
    def reseq_human_files(self):
        return [("hq_reseq_to_hg_bam", self.hq_reseq_to_hg_bam),
                ("flnc_reseq_to_hg_bam", self.flnc_reseq_to_hg_bam),
                ("hq_reseq_to_hg_selected_transcripts_csv", self.hq_reseq_to_hg_selected_transcripts_csv),
                ("flnc_reseq_to_hg_selected_transcripts_csv", self.flnc_reseq_to_hg_selected_transcripts_csv)]

    @property
    def sirv_files(self):
        """sirv related files."""
        return [
            ("collapsed_to_sirv_rep_fastq", self.collapsed_to_sirv_rep_fq),
            ("collapsed_to_sirv_rep_bam", self.collapsed_to_sirv_rep_bam),
            ("collapsed_to_sirv_rep_readlength_csv", self.collapsed_to_sirv_rep_readlength_csv),
            ("chained_ids_txt", self.chained_ids_txt)
        ]

    def ln_gencode_gtf(self, gencode_gtf):
        """Make link of gencode_gtf"""
        log.info("Making soft link of gencode annotation gtf")
        lnabs(gencode_gtf, self.gencode_gtf)

    def ln_sirv_truth_dir(self, sirv_truth_dir):
        """Make a link of sirv ground truth dir"""
        log.info("Making soft link of sirv ground truth dir")
        lnabs(sirv_truth_dir, self.sirv_truth_dir)

    def make_all_files_from_SMRTLink_job(self,  make_readlength=True):
        """Make all data files from a smrtlink job, including
        * consolidate flnc and nfl xml files to fasta
        * collect ccs, classify, cluster reports from sl job and make validation_report_csv
        * make read length csv files
        """
        mkdir(self.root_dir)
        mkdir(self.fasta_dir)
        mkdir(self.csv_dir)
        mkdir(self.chain_sample_dir)
        mkdir(self.reseq_to_sirv_dir)
        mkdir(self.reseq_to_hg_dir)

        smrtlink_job_dir = self.smrtlink_job_dir
        consolidate_xml(src=self.sl_job.flnc_bam, dst=self.isoseq_flnc_fa)

        # symlink smrtlink_job_dir and files to validation dir
        self.ln_files_from_SMRTLink_job()
        self.make_readlength_csvs()

    def make_readlength_csvs(self):
        """Make all read length csv files."""
        log.info("make all readlength csv files.")
        z = [
            (self.ccs_fa, self.ccs_readlength_csv),
            (self.isoseq_flnc_fa, self.flnc_readlength_csv),
            (self.hq_isoforms_fa, self.hq_readlength_csv),
            (self.lq_isoforms_fa, self.lq_readlength_csv),
            (self.consensus_isoforms_fa, self.consensus_isoforms_readlength_csv)
        ]
        for fasta_fn, csv_fn in z:
            make_readlength_csv(fasta_fn=fasta_fn, csv_fn=csv_fn)

    def make_readlength_csv_for_sirv_isoforms(self):
        """Make read length csv for representative isoforms collapsing to sirv"""
        make_readlength_csv(fasta_fn=self.collapsed_to_sirv_rep_fq,
                            csv_fn=self.collapsed_to_sirv_rep_readlength_csv)

    def make_readlength_csv_for_hg_isoforms(self):
        """Make read length csv for representative isoforms collapsing to hg"""
        make_readlength_csv(fasta_fn=self.collapsed_to_hg_rep_fq,
                            csv_fn=self.collapsed_to_hg_rep_readlength_csv)

    def ln_files_from_SMRTLink_job(self):
        """Make soft links to existing isoseq output fasta|fastq files."""
        log.info("make soft links from smrtlink job")
        # Make a link of smrtlink dir
        lnabs(self.smrtlink_job_dir, op.join(self.root_dir, op.basename(self.sl_job.root_dir)))

        # Make a link of consensus isoforms fa, hq|lq isoforms fa|fq, isoseq_flnc.fasta
        lnabs(src=self.sl_job.hq_isoforms_fa, dst=self.hq_isoforms_fa)
        lnabs(src=self.sl_job.hq_isoforms_fq, dst=self.hq_isoforms_fq)
        lnabs(src=self.sl_job.lq_isoforms_fa, dst=self.lq_isoforms_fa)
        lnabs(src=self.sl_job.lq_isoforms_fq, dst=self.lq_isoforms_fq)

        if op.exists(self.sl_job.ccs_xml):
            consolidate_xml(self.sl_job.ccs_xml, self.ccs_fa)
        else:
            raise IOError("Could neither find {}".format(self.sl_job.ccs_xml))

        self.sl_job.export_unpolished_fa(unpolished_fa=self.consensus_isoforms_fa)

    def reseq_to_human(self, target_fa, selected_transcripts):
        """Resequence FLNC and HQ isoforms to human transcripts"""
        log.info('Locals: %s' % locals())
        _files = [(self.hq_isoforms_fa, target_fa, self.hq_reseq_to_hg_bam, self.hq_reseq_to_hg_selected_transcripts_csv),
                  (self.isoseq_flnc_fa, target_fa, self.flnc_reseq_to_hg_bam, self.flnc_reseq_to_hg_selected_transcripts_csv)]
        for query_fa, target_fa, out_bam, out_csv in _files:
            reseq(query_fa, target_fa, out_bam, 16)
            with open(out_csv, 'w') as writer:
                writer.write(coverage2str(coverage_d=subset_dict(
                    d=bam2coverage(out_bam), selected_keys=selected_transcripts)))

    def make_readme_txt(self, args):
        """Write args and data files to README file."""
        with open(self.readme_txt, 'w') as writer:
            log.info("args=%s\n", args)
            writer.write("# Created by pbtranscript-internal-validation.ValidationRunner.make_readme_txt()\n")
            writer.write("args=%s\n\n" % args)

            files = self.common_files + self.collapse_human_files + self.reseq_human_files + self.sirv_files
            for desc, fn in files:
                if op.exists(fn):
                    writer.write("%s=%s\n" % (desc, fn))
