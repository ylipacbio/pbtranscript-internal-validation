#!/bin/bash
source bamboo_setup.sh
set -vex

# debug bamboo
ls -larth
ls -larth $THISDIR/..

make clean

WHEELHOUSE=/home/cdunn/wheelhouse/gcc-6
ls -larth ${WHEELHOUSE}

# Build with dependencies (fairly fast)
pip install --user --find-links=${WHEELHOUSE} --no-index pytest pytest-cov pylint cython


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

if [ -e $THISDIR/../pbcommand ] ; then
    pushd ../pbcommand
    pwd
    ls
    pip install -v --user --find-links=${WHEELHOUSE} --no-index --no-deps --edit  .
    popd
else
    pip install -v --user --find-links=${WHEELHOUSE} --no-index pbcommand
fi

if [ -e $THISDIR/../pbcoretools ] ; then
    pushd ../pbcoretools
    pip install -v --user --find-links=${WHEELHOUSE} --no-index --no-deps --edit  .
    popd
else
    pip install -v --user --find-links=${WHEELHOUSE} --no-index pbcoretools
fi


if [ -e $THISDIR/../pbtranscript ] ; then
    pushd ../pbtranscript
    pip install -v --user --find-links=${WHEELHOUSE} --no-index --no-deps --edit .
    popd
else
    pip install -v --user --find-links=${WHEELHOUSE} --no-index pbtranscript
fi

if [ -e $THISDIR/../pbtranscript2 ] ; then
    pushd ../pbtranscript2
    pip install -v --user --find-links=${WHEELHOUSE} --no-index --no-deps --edit .
    popd
else
    pip install -v --user --find-links=${WHEELHOUSE} --no-index pbtranscript2
fi

python -c "import pbcore"
python -c "import pbcoretools"
python -c "import pbcommand"
python -c "import pbtranscript"
python -c "import pbtranscript2"

pip install -v --user --edit . 

# Test.
make pylint
time make test

pwd
ls -larth
