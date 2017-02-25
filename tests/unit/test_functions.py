#!/usr/bin/env python

import unittest
import os.path as op
import filecmp
from pbcore.util.Process import backticks
from pbtranscript.Utils import mkdir, rmpath
from pbtranscript_internal_validation.validate_smrtlink_isoseq_rc0 import get_gtf_chromosomes, filter_sam_by_targets

from test_setpath import smrtlink_isoseq_jobs, OUT_DIR


def test_get_gtf_chrsomsomes():
    gtf_chrs = get_gtf_chromosomes('/pbi/dept/secondary/siv/testdata/isoseq/gencode/gencode.v25.annotation.gtf')
    assert gtf_chrs == set(["chr%s" % x for x in range(1, 23) + ['X', 'Y', 'M']])
    return gtf_chrs


def test_filter_sam_by_targets():
    d = op.join(OUT_DIR, 'test_filter_sam_by_targets')
    rmpath(d)
    mkdir(d)
    in_sam = op.join(d, 'in.sam')
    out_sam = op.join(d, 'out.sam')
    filtered_sam = op.join(d, 'filted.sam')
    print in_sam
    print out_sam
    print filtered_sam

    comment = '#random\n'
    aln1 = 'read1\tsomething\tchr1\n'
    aln2 = 'read2\tsomething\tutg_scaffold\n'

    with open(in_sam, 'w') as writer:
        writer.write(comment + aln1 + aln2)

    filter_sam_by_targets(in_sam=in_sam, targets=test_get_gtf_chrsomsomes(),
                          out_sam=out_sam, filtered_sam=filtered_sam)
    op.exists(out_sam)
    op.exists(filtered_sam)
    for l in open(out_sam, 'r'):
        print l

    assert [l for l in open(out_sam, 'r')] == [comment, aln1]
    assert [l for l in open(filtered_sam, 'r')] == [comment, aln2]
