from __future__ import print_function
from setuptools import setup, find_packages
from distutils.extension import Extension
import os.path
import sys
import numpy

__author__ = "etseng|yli@pacificbiosciences.com"
version = "3.0.0"

if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__


def _get_local_file(file_name):
    return os.path.join(os.path.dirname(__file__), file_name)


def _get_requirements(file_name):
    with open(file_name, 'r') as f:
        reqs = [line for line in f if not line.startswith("#")]
    return reqs


def _get_local_requirements(file_name):
    return _get_requirements(_get_local_file(file_name))


def run_cmd(cmd):
    import subprocess
    p = subprocess.Popen(cmd,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         shell=True)
    return p.stdin, p.stdout

def check_program(program):
    import re
    try:
        stdin, stdout = run_cmd('/usr/bin/which "%s"||echo no' % program)
        lines = stdout.readlines()
        match = filter(lambda x: re.search(program, x), lines)
        if match:
            return True
        else:
            return False
    except (IndexError, ValueError, AttributeError):
        print('%s is required, please install %s' % program, file=sys.stderr)
        return False


def exit_if_not_installed(program):
    if (not check_program(program=program)):
        print('Unable to install - %s must be installed and in the PATH variable' % program, file=sys.stderr)
        sys.exit(1)


if ("install" in sys.argv):
    exit_if_not_installed("pbtranscript")

setup(
    name = 'pbtranscript_internal_validation',
    version=version,
    author='Pacific Biosciences',
    author_email='devnet@pacificbiosciences.com',
    license='LICENSE.txt',
    ext_modules = [],
    include_dirs=[numpy.get_include()],
    scripts=[], #'pbtranscript_internal_validation/validate_smrtlink_isoseq_rc0.py'],
    entry_points={'console_scripts': [
        'validate_smrtlink_isoseq_rc0.py = pbtranscript_internal_validation.validate_smrtlink_isoseq_rc0:main',
        'validate_smrtlink_isoseq2_rc0.py = pbtranscript_internal_validation.validate_smrtlink_isoseq2_rc0:main',
    ]},
    package_dir={'pbtranscript_internal_validation': 'pbtranscript_internal_validation'},
    #package_data={'pbtranscript_internal_validation': []},
    packages=find_packages(),
    zip_safe=False,
    install_requires=_get_local_requirements("REQUIREMENTS.txt")
)
