from desc.sims.GCRCatSimInterface import PhoSimDESCQA
from desc.sims.GCRCatSimInterface import CompoundDESCQAInstanceCatalog
from desc.sims.GCRCatSimInterface import bulgeDESCQAObject_protoDC2
from desc.sims.GCRCatSimInterface import diskDESCQAObject_protoDC2
from desc.sims.GCRCatSimInterface import CompoundDESCQAObject
from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap
from lsst.sims.utils import arcsecFromRadians
from lsst.sims.utils import ObservationMetaData
import numpy as np
import os

class _testDESCQAObj(object):
    yaml_file_name = 'protoDC2'
    field_ra = 23.0
    field_dec = -22.3


# use different DESCQAObject classes with different
# _cat_cache_suffixes for baseline and test to force
# the CompoundDESCQAInstanceCatalog to load the catalog
# and apply the rotation itself

class bulgeDESCQAObject_baseline(_testDESCQAObj, bulgeDESCQAObject_protoDC2):
    _cat_cache_suffix = 'rot_baseline'

class diskDESCQAObject_baseline(_testDESCQAObj, diskDESCQAObject_protoDC2):
    _cat_cache_suffix = 'rot_baseline'

class bulgeDESCQAObject_test(_testDESCQAObj, bulgeDESCQAObject_protoDC2):
    objid = 'bulge_descqa'
    _cat_cache_suffix = 'rot_test'

class diskDESCQAObject_test(_testDESCQAObj, diskDESCQAObject_protoDC2):
    objid = 'disk_descqa'
    _cat_cache_suffix = 'rot_test'

class PhoSimDESCQABulge(PhoSimDESCQA):
    cannot_be_null = ['hasBulge']

class PhoSimDESCQADisk(PhoSimDESCQA):
    cannot_be_null = ['hasDisk']


if __name__ == "__main__":

    obs = ObservationMetaData(pointingRA=22.7, pointingDec=-22.7,
                              boundType='circle', boundLength=0.01,
                              mjd=59981.2, bandpassName='g',
                              rotSkyPos=18.6)

    out_dir = 'compound_verification_dir'
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    else:
        if not os.path.isdir(out_dir):
            raise RuntimeError('%s is not a dir' % out_dir)

    # first, write out the two InstanceCatalogs separately
    bulge_db = bulgeDESCQAObject_baseline()
    bulge_cat = PhoSimDESCQABulge(bulge_db, obs_metadata=obs)
    bulge_cat_name = os.path.join(out_dir, 'bulge_baseline.txt')
    bulge_cat.phoSimHeaderMap = DefaultPhoSimHeaderMap
    bulge_cat.write_catalog(bulge_cat_name, chunk_size=10000)

    disk_db = diskDESCQAObject_baseline()
    disk_cat = PhoSimDESCQADisk(disk_db, obs_metadata=obs)
    disk_cat_name = os.path.join(out_dir, 'disk_baseline.txt')
    disk_cat.phoSimHeaderMap = DefaultPhoSimHeaderMap
    disk_cat.write_catalog(disk_cat_name, chunk_size=10000)

    baseline_lines = set()
    with open(bulge_cat_name, 'r') as in_file:
        for line in in_file:
            baseline_lines.add(line)

    assert len(baseline_lines) > 100
    i0 = len(baseline_lines)

    with open(disk_cat_name, 'r') as in_file:
        for line in in_file:
            baseline_lines.add(line)

    assert len(baseline_lines) > (i0+100)

    # now, write the same catalogs with the CompoundDESCQAInstanceCatalog
    cat = CompoundDESCQAInstanceCatalog([PhoSimDESCQABulge,
                                         PhoSimDESCQADisk],
                                        [bulgeDESCQAObject_test,
                                         diskDESCQAObject_test],
                                        obs_metadata=obs,
                                        compoundDBclass=CompoundDESCQAObject)

    cat.phoSimHeaderMap = DefaultPhoSimHeaderMap

    compound_name = os.path.join(out_dir, 'compound_cat.txt')
    cat.write_catalog(compound_name, chunk_size=10000)

    n_rows = 0
    with open(compound_name, 'r') as in_file:
        for line in in_file:
            try:
                assert line in baseline_lines
            except AssertionErrro:
                print('\n\n%s\n\nmissing\n' % line)
                raise
            n_rows += 1

    # because baseline_lines is a set and sets do not allow
    # duplicate entries, the header lines, which should be
    # identical in both the bulge and the disk catalogs,
    # will only appear once in baseline_lines
    assert n_rows == len(baseline_lines)
