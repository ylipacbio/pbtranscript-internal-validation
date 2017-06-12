#!/usr/bin/env python

import unittest
import os.path as op
import filecmp
from pbcore.util.Process import backticks
from pbtranscript.Utils import mkdir, rmpath
from pbtranscript_internal_validation.Utils import *

from test_setpath import smrtlink_isoseq_jobs, OUT_DIR


m4_str = \
"""qName tName score percentSimilarity qStrand qStart qEnd qLength tStrand tStart tEnd tLength mapQV
i0_HQ_sample60b8a7|c50/f7p0/839 MT-ATP6-201|MT -3352 99.5569 0 162 839 839 0 0 677 681 254
i0_HQ_sample60b8a7|c42/f6p0/797 PFN1-001|PFN1 -3960 99.3766 0 0 797 797 0 484 1286 1291 254
i0_HQ_sample60b8a7|c42/f6p0/797 PFN1-001|PFN1 -3960 99.3766 0 0 797 797 0 484 1286 1291 254
i0_HQ_sample60b8a7|c42/f6p0/797 PFN1-001|PFN1 -3960 99.3766 0 0 797 797 0 484 1286 1291 254
i0_HQ_sample60b8a7|c52/f7p0/707 MT-CO2-201|MT -3398 99.7076 0 0 684 707 0 0 684 684 254
i0_HQ_sample60b8a7|c16/f3p0/882 IGLC2-001|IGLC2 -2228 99.3363 0 431 882 882 0 0 452 462 0
i0_HQ_sample60b8a7|c44/f3p0/869 FTL-001|FTL -4299 99.4253 0 0 869 869 0 9 876 878 254"""


def assert_dict_equal(d, expected):
    assert sorted(expected.keys()) == sorted(d.keys())
    for key in expected.keys():
        assert d[key] == expected[key]

def test_m42coverage():
    out_m4 = op.join(OUT_DIR, 'test_m42coverage.m4')
    with open(out_m4, 'w') as writer:
        writer.write(m4_str)
    d = m42coverage(out_m4)
    expected = {'MT-ATP6-201|MT':1, 'PFN1-001|PFN1':3, 'MT-CO2-201|MT':1, 'IGLC2-001|IGLC2':1, 'FTL-001|FTL':1}
    assert_dict_equal(d, expected)

def test_subset_dict():
    in_d = {'MT-ATP6-201|MT':1, 'PFN1-001|PFN1':3, 'MT-CO2-201|MT':1, 'IGLC2-001|IGLC2':1, 'FTL-001|FTL':1}
    expected = {'PFN1-001|PFN1':3, 'MT-CO2-201|MT':1}
    assert_dict_equal(subset_dict(in_d, expected.keys()), expected)

def coverage2str():
    in_d = {'MT-ATP6-201|MT':1, 'PFN1-001|PFN1':3}
    assert covearge2str(in_d) == 'transcript\tcoverage\nMT-ATP6-201|MT\t1\nPFN1-001|PFN1\t3\n'
    #coverage2str(coverage_d=subset_dict(d=m42coverage(out_m4), selected_keys=selected_transcripts)
