#!/bin/bash
source scl_source enable devtoolset-8
source loadLSST.bash
setup -t sims_w_2019_42 lsst_sims
setup -t DC2production throughputs
setup -t DC2production sims_skybrightness_data
pip install nose
pip install pylint
pip install gcr
git clone https://github.com/LSSTDESC/gcr-catalogs.git
export PYTHONPATH=`pwd`/gcr-catalogs:${PYTHONPATH}
eups declare sims_GCRCatSimInterface -r ${TRAVIS_BUILD_DIR} -t current
setup sims_GCRCatSimInterface
cd ${TRAVIS_BUILD_DIR}
scons
nosetests -s --with-coverage --cover-package=desc.sims.GCRCatSimInterface
