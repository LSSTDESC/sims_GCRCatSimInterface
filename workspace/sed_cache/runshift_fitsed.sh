#!/bin/bash
source /opt/lsst/software/stack/loadLSST.bash
setup lsst_distrib
setup lsst_sims
which python
eups list lsst_distrib
