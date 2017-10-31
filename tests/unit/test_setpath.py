#!/usr/bin/python
from os import path
import ConfigParser

"""Define test data path for pbtranscript-internal-validation."""

THIS_DIR = path.dirname(path.abspath(__file__))
ROOT_DIR = path.dirname(THIS_DIR)
NOSE_CFG = path.join(THIS_DIR, "nose.cfg")

def _get_smrtlink_isoseq_jobs():
    """Get smrtlink isoseq jobs for unittest"""
    nosecfg = ConfigParser.SafeConfigParser()
    nosecfg.readfp(open(NOSE_CFG), 'r')
    if nosecfg.has_section('smrtlink_isoseq_jobs'):
        nti = path.abspath(nosecfg.get('smrtlink_isoseq_jobs', 'NTI'))
        rc0 = path.abspath(nosecfg.get('smrtlink_isoseq_jobs', 'RC0'))
        tiny = path.abspath(nosecfg.get('smrtlink_isoseq_jobs', 'TINY')) # tiny rc0
        return {'NTI': nti, 'RC0': rc0, 'TINY': tiny}
    else:
        msg = "Unable to find section [smrtlink_isoseq_jobs] option [NTI] in {f}".format(f=NOSE_CFG)
        raise KeyError(msg)

OUT_DIR = path.join(ROOT_DIR, "out")
DATA_DIR = path.join(ROOT_DIR, "data")
STD_DIR = path.join(ROOT_DIR, "stdout")
smrtlink_isoseq_jobs = _get_smrtlink_isoseq_jobs()
