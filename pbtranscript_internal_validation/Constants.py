#!/usr/bin/env python

"""Define constants of data resources to use for validation."""
import os.path as op


class ValidationConstants(object):

    """ Constants of validation resources and parameter constants"""
    NPROC = 16
    SIRV_DATA_ROOT_DIR = "/pbi/dept/secondary/siv/testdata/isoseq/lexigoen-ground-truth"

    GENCODE_ROOT_DIR = "/pbi/dept/secondary/siv/testdata/isoseq/gencode"
    GENCODE_GTF = op.join(GENCODE_ROOT_DIR, "gencode.v25.annotation.gtf")
    HUMAN_TRANSCRIPTS_FA = op.join(GENCODE_ROOT_DIR, "gencode.v25.transcripts.fa")
    HUMAN_SIRV_REFERENCE = "/pbi/dept/secondary/siv/testdata/isoseq/references/hg38_and_sirv/hg38_and_sirv/referenceset.xml"
    HUMAN_REFERENCE = "/pbi/dept/secondary/siv/references/human_GRCh38_no_alt_analysis_set/sequence/human_GRCh38_no_alt_analysis_set.fasta"
    HUMAN_12_TRANSCRIPTS_FA = op.join(SIRV_DATA_ROOT_DIR, 'human', 'RC0_selected_12_human_transcripts.fasta')
    HUMAN_16_TRANSCRIPTS_FA = op.join(SIRV_DATA_ROOT_DIR, 'human', 'RC0_selected_16_human_transcripts.fasta')

    SIRV_NAME = 'SIRV'
    SIRV_TRUTH_DIR = op.join(SIRV_DATA_ROOT_DIR, 'validation')
    SIRV_TRANSCRIPTS_FA = op.join(SIRV_DATA_ROOT_DIR, 'isoforms', 'SIRV_C_150601a.truth.fasta')
    SIRV_REFERENCE = op.join(SIRV_DATA_ROOT_DIR, 'reference', 'SIRV_150601a.fasta')
