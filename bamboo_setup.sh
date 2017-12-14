#!/bin/bash
set -vex
type module >& /dev/null || . /mnt/software/Modules/current/init/bash

# This must come first so that everything else we add later
# will take precedence.
module purge
module load smrttools/incremental
which blasr
which pbindex
which pbdagcon

module load gcc
module load zlib
module load ccache
module load python/2

# For isolation:
export PYTHONUSERBASE=$(pwd)/LOCAL
mkdir -p LOCAL/
export PATH=${PYTHONUSERBASE}/bin:${PATH}
which python
which pip

THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo THISDIR=$THISDIR
cd $THISDIR
# Everything else happens in THISDIR.
