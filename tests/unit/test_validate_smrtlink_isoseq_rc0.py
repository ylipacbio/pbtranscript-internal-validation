#!/usr/bin/env python

import unittest
import os.path as op
import filecmp
from pbcore.util.Process import backticks
from pbtranscript_internal_validation import ValidationFiles
from pbtranscript_internal_validation.validate_smrtlink_isoseq_rc0 import main

from test_setpath import smrtlink_isoseq_jobs, OUT_DIR

def test_cmd_call():
    """Test calling validate_smrtlink_isoseq_rc0.py from command line"""
    print smrtlink_isoseq_jobs
    nti_dir = smrtlink_isoseq_jobs['NTI']
    eval_dir = OUT_DIR

    #cmd = "validate_smrtlink_isoseq_rc0.py %s %s" % (nti_dir, eval_dir)
    #print cmd

    main(args=[nti_dir, eval_dir, '--sirv_only'])

    #m, c, e = backticks(cmd)
    #assert c == 0

    #for desc, fn in ValidationFiles(eval_dir).all_files:
    #    assert op.exists(fn)

