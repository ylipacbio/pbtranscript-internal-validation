#!/usr/bin/env python

"""Define constants of data resources to use for validation."""
import os.path as op


__all__ = [ "ValidationConstants" ]


class ValidationConstants(object):

    """ Constants of validation resources and parameter constants"""
    GENCODE_ROOT_DIR = "/pbi/dept/secondary/siv/testdata/isoseq/gencode"
    GENCODE_GTF = op.join(GENCODE_ROOT_DIR, "gencode.v25.annotation.gtf")
    HUMAN_TRANSCRIPTS_FA = op.join(GENCODE_ROOT_DIR, "gencode.v25.transcripts.fa")

    GMAP_DB = "/pbi/dept/secondary/siv/testdata/isoseq/gmap_db/"
    HUMAN_GMAP_NAME = "hg38"
    SIRV_GMAP_NAME = "SIRV"
    GMAP_NPROC = 12
    SIRV_NAME = "SIRV"
    SIRV_DATA_ROOT_DIR = "/pbi/dept/secondary/siv/testdata/isoseq/lexigoen-ground-truth"
    SIRV_TRUTH_DIR = op.join(SIRV_DATA_ROOT_DIR, 'validation')
    SIRV_TRANSCRIPTS_FA = op.join(SIRV_DATA_ROOT_DIR, 'isoforms', 'SIRV_C_150601a.truth.fasta')
    SIRV_GENOME_FA = op.join(SIRV_DATA_ROOT_DIR, 'reference', 'SIRV_150601a.fasta')
