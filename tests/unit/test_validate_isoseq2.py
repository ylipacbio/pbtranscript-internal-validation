#!/usr/bin/env python
"""
Validate an IsoSeq SMRTLink job of RC0 sample.
"""
import unittest
import os.path as op
import filecmp
from pbcore.util.Process import backticks
from pbtranscript.Utils import mkdir, rmpath
from pbtranscript_internal_validation import ValidationFiles, ValidationRunner
from pbtranscript_internal_validation.validate_smrtlink_isoseq2_rc0 import main

from test_setpath import smrtlink_isoseq_jobs, OUT_DIR


def test_validate_smrtlink_isoseq_rc0():
    """Test calling validate_smrtlink_isoseq_rc0.py from command line"""
    print smrtlink_isoseq_jobs
    rc0_dir = smrtlink_isoseq_jobs['TINY']
    eval_dir = op.join(OUT_DIR, 'test_validate_rc0')
    rmpath(eval_dir)
    mkdir(eval_dir)

    # identical to "validate_smrtlink_isoseq2_rc0.py %s %s --collapse_to_human --reseq_to_human --make_readlength" % (rc0_dir, eval_dir)
    #main(args=[rc0_dir, eval_dir, '--collapse_to_human', '--reseq_to_human', '--make_readlength'])

    # identical to "validate_smrtlink_isoseq2_rc0.py %s %s" % (rc0_dir, eval_dir)
    main(args=[rc0_dir, eval_dir, '--collapse_to_human', '--reseq_to_human', '--make_readlength'])

    runner = ValidationRunner(eval_dir, rc0_dir)
    for desc, fn in runner.all_files:
        if not op.exists(fn):
            print 'File %s does not exist' % fn
            assert 'File %s does not exist' % fn == False
