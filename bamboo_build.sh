#!/bin/bash
source bamboo_setup.sh
set -vex

# debug bamboo
ls -larth
ls -larth $THISDIR/..

# If .pyc exists from previous run, and if files have moved, tests can fail. #TODO(CD): drop this
#find pbtranscript2 -name '*.pyc' | xargs rm -f


# Reinvent the wheel (7s -- turn back on someday maybe)
#python setup.py bdist_wheel
# (goes into ./dist/)
# (not needed yet, but a useful quick check)

WHEELHOUSE=/home/cdunn/wheelhouse/gcc-6
ls -larth ${WHEELHOUSE}


pip install -v --user --find-links=${WHEELHOUSE} --no-index pytest
make test-fast

# Build with dependencies (fairly fast)
pip install --user --find-links=${WHEELHOUSE} --no-index pytest pytest-cov pylint cython

#pip --no-cache-dir install ... ?

ZLIB_CFLAGS=$(pkg-config zlib --cflags)
ZLIB_LIBS=$(pkg-config zlib --libs)
HTSLIB_CONFIGURE_OPTIONS="--disable-libcurl --disable-bz2 --disable-lzma" \
  CFLAGS="-D_GNU_SOURCE ${ZLIB_CFLAGS}" \
  LDFLAGS="${ZLIB_LIBS}" \
  pip install -v --user --find-links=${WHEELHOUSE} --no-index 'pysam'
#  pip install -v --user 'pysam==0.9.1.4'
python -c 'import pysam.version; print pysam.version.__version__'
python -c 'import pysam.cfaidx; print pysam.cfaidx' || python -c 'import pysam.libcfaidx; print pysam.libcfaidx'

python -c 'import numpy; print numpy.__version__' ||  pip install -v --user 'numpy==1.9'
pip install -v --user 'networkx==1.1'
python -c "import networkx; print networkx.__version__;"
python -c "import networkx; a=networkx.Graph(); a.add_edge('a', 'b'); a.degree().items();"

if [ -e $THISDIR/../pbcore ] ; then
    pushd ../pbcore
    pip install -v --user --find-links=${WHEELHOUSE} --no-index --edit .
    popd
else
    pip install -v --user --find-links=${WHEELHOUSE} --no-index pbcore
fi

# Drop --no-index, in case we are missing something.
pip install -v --user --find-links=${WHEELHOUSE} --edit .

# Test.
export MY_TEST_FLAGS="-v -s --durations=0 --cov=pbtranscript2 --cov-report=term-missing --cov-report=xml:coverage.xml --cov-branch"
time make pytest
sed -i -e 's@filename="@filename="./pbtranscript2/@g' coverage.xml

time make pylint


pwd
ls -larth
