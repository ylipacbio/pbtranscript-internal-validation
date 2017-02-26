#!/usr/bin/env python
"""
Validate an IsoSeq SMRTLink job of human caner NTI sample.
NOTE: This validation is slow
"""

import unittest
import os.path as op
import filecmp
from pbcore.util.Process import backticks
from pbtranscript.Utils import mkdir, rmpath
from pbtranscript_internal_validation import ValidationFiles, ValidationRunner
from pbtranscript_internal_validation.validate_smrtlink_isoseq_rc0 import main

from test_setpath import smrtlink_isoseq_jobs, OUT_DIR


def test_validate_nti():
    """Test calling validate_smrtlink_isoseq_rc0.py from command line"""
    print smrtlink_isoseq_jobs
    nti_dir = smrtlink_isoseq_jobs['NTI']
    eval_dir = op.join(OUT_DIR, 'test_validate_nti')
    rmpath(eval_dir)
    mkdir(eval_dir)

    # identical to "validate_smrtlink_isoseq_rc0.py %s %s" % (nti_dir, eval_dir)
    main(args=[nti_dir, eval_dir, '--human_only'])

    runner = ValidationRunner(eval_dir, nti_dir)
    for desc, fn in runner.common_files + runner.human_files:
        if not op.exists(fn):
            print 'File %s does not exist' % fn
            assert 'File %s does not exist' % fn == False
