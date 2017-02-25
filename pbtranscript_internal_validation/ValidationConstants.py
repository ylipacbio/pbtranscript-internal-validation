#!/usr/bin/env python

"""Define constants of data resources to use for validation."""


__all__ = [ "ValidationConstants" ]


class ValidationConstants(object):

    """ Constants of validation resources and parameter constants"""

    GENCODE_GTF = "/pbi/dept/secondary/siv/testdata/isoseq/gencode/gencode.v25.annotation.gtf"
    GMAP_DB = "/pbi/dept/secondary/siv/testdata/isoseq/gmap_db/"
    HUMAN_GMAP_NAME = "hg38"
    SIRV_GMAP_NAME = "SIRV"
    GMAP_NPROC = 12
    SIRV_NAME = "SIRV"
    SIRV_TRUTH_DIR = "/pbi/dept/secondary/siv/testdata/isoseq/lexigoen-ground-truth/validation"
