"""
Code to write instance catalogs for phosim and imsim using the gcr-catalogs
galaxy catalogs.
"""
from __future__ import with_statement
import os
import copy
import subprocess
from collections import namedtuple
import gzip
import numpy as np
import h5py

from lsst.sims.catalogs.definitions import parallelCatalogWriter
from lsst.sims.catalogs.decorators import cached
from lsst.sims.catUtils.baseCatalogModels import StarObj
from lsst.sims.catUtils.exampleCatalogDefinitions import \
    PhoSimCatalogPoint, DefaultPhoSimHeaderMap
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.utils import arcsecFromRadians, _getRotSkyPos
from . import PhoSimDESCQA, bulgeDESCQAObject, diskDESCQAObject, knotsDESCQAObject

__all__ = ['InstanceCatalogWriter', 'make_instcat_header', 'get_obs_md']

class InstanceCatalogWriter(object):
    """
    Class to write instance catalogs for PhoSim and imSim using
    galaxies accessed via the gcr-catalogs interface.
    """
    def __init__(self, opsimdb, descqa_catalog, dither=True,
                 min_mag=10, minsource=100, proper_motion=False,
                 imsim_catalog=False):
        """
        Parameters
        ----------
        obsimdb: str
            OpSim db filename.
        descqa_catalog: str
            Name of the DESCQA galaxy catalog.
        dither: bool [True]
            Flag to enable the dithering included in the opsim db file.
        min_mag: float [10]
            Minimum value of the star magnitude at 500nm to include.
        minsource: int [100]
            Minimum number of objects for phosim.py to simulate a chip.
        proper_motion: bool [True]
            Flag to enable application of proper motion to stars.
        imsim_catalog: bool [False]
            Flag to write an imsim-style object catalog.
        """
        if not os.path.exists(opsimdb):
            raise RuntimeError('%s does not exist' % opsimdb)

        self.descqa_catalog = descqa_catalog
        self.dither = dither
        self.min_mag = min_mag
        self.minsource = minsource
        self.proper_motion = proper_motion
        self.imsim_catalog = imsim_catalog

        self.obs_gen = ObservationMetaDataGenerator(database=opsimdb,
                                                    driver='sqlite')

        self.star_db = StarObj(database='LSSTCATSIM',
                               host='fatboy.phys.washington.edu',
                               port=1433, driver='mssql+pymssql')

        self.instcats = get_instance_catalogs(imsim_catalog)

    def write_catalog(self, obsHistID, out_dir='.', fov=2):
        """
        Write the instance catalog for the specified obsHistID.

        Parameters
        ----------
        obsHistID: int
            ID of the desired visit.
        out_dir: str ['.']
            Output directory.  It will be created if it doesn't already exist.
        fov: float [2.]
            Field-of-view angular radius in degrees.  2 degrees will cover
            the LSST focal plane.
        """
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        obs_md = get_obs_md(self.obs_gen, obsHistID, fov, dither=self.dither)

        cat_name = 'phosim_cat_%d.txt' % obsHistID
        star_name = 'star_cat_%d.txt' % obsHistID
        bright_star_name = 'bright_stars_%d.txt' % obsHistID
        gal_name = 'gal_cat_%d.txt' % obsHistID
        knots_name = 'knots_cat_%d.txt' % obsHistID
        #agn_name = 'agn_cat_%d.txt' % obshistid

        make_instcat_header(self.star_db, obs_md, cat_name,
                            imsim_catalog=self.imsim_catalog,
                            object_catalogs=(star_name, gal_name))

        star_cat = self.instcats.StarInstCat(self.star_db, obs_metadata=obs_md,
                                             cannot_be_null=['inProtoDc2'])
        star_cat.min_mag = self.min_mag
        star_cat.disable_proper_motion = not self.proper_motion

        bright_cat \
            = self.instcats.BrightStarInstCat(self.star_db, obs_metadata=obs_md,
                                              cannot_be_null=['isBright'])
        bright_cat.min_mag = self.min_mag

        cat_dict = {os.path.join(out_dir, star_name): star_cat,
                    os.path.join(out_dir, bright_star_name): bright_cat}
        parallelCatalogWriter(cat_dict, chunk_size=100000, write_header=False)
        
        # TODO: Find a better way of checking for catalog type
        if 'knots' in self.descqa_catalog:
            cat = self.instcats.DESCQACat(knotsDESCQAObject(self.descqa_catalog),
                                          obs_metadata=obs_md,
                                          cannot_be_null=['hasKnots'])
            cat.write_catalog(os.path.join(out_dir, knots_name), chunk_size=100000,
                              write_header=False)
        else:
            # Creating empty knots component
            subprocess.check_call('cd %(out_dir)s; touch %(knots_name)s' % locals(), shell=True)
            
        cat = self.instcats.DESCQACat(bulgeDESCQAObject(self.descqa_catalog),
                                      obs_metadata=obs_md,
                                      cannot_be_null=['hasBulge'])
        cat.write_catalog(os.path.join(out_dir, gal_name), chunk_size=100000,
                          write_header=False)

        cat = self.instcats.DESCQACat(diskDESCQAObject(self.descqa_catalog),
                                      obs_metadata=obs_md,
                                      cannot_be_null=['hasDisk'])
        cat.write_catalog(os.path.join(out_dir, gal_name), chunk_size=100000,
                          write_mode='a', write_header=False)

        if self.imsim_catalog:
            
            imsim_cat = 'imsim_cat_%i.txt' % obsHistID
            command = 'cd %(out_dir)s; cat %(cat_name)s %(star_name)s %(gal_name)s %(knots_name) > %(imsim_cat)s' % locals()
            subprocess.check_call(command, shell=True)

        # gzip the object files.
        for orig_name in (star_name, gal_name):
            full_name = os.path.join(out_dir, orig_name)
            with open(full_name, 'rb') as input_file:
                with gzip.open(full_name+'.gz', 'wb') as output_file:
                    output_file.writelines(input_file)
            os.unlink(full_name)


def make_instcat_header(star_db, obs_md, outfile, object_catalogs=(),
                        imsim_catalog=False, nsnap=1, vistime=30.,
                        minsource=100):
    """
    Write the header part of an instance catalog.

    Parameters
    ----------
    star_db: lsst.sims.catUtils.baseCatalogModels.StarObj
        InstanceCatalog object for stars.  Connects to the UW fatboy db server.
    obs_md: lsst.sims.utils.ObservationMetaData
        Observation metadata object.
    object_catalogs: sequence [()]
        Object catalog names to include in base phosim instance catalog.
        Defaults to an empty tuple.
    imsim_catalog: bool [False]
        Flag to write an imSim-style object catalog.
    nsnap: int [1]
        Number of snaps per visit.
    vistime: float [30.]
        Visit time in seconds.
    minsource: int [100]
        Minimum number of objects for phosim.py to simulate a chip.

    Returns
    -------
    lsst.sims.catUtils.exampleCatalogDefinitions.PhoSimCatalogPoint object
    """
    cat = PhoSimCatalogPoint(star_db, obs_metadata=obs_md)
    cat.phoSimHeaderMap = copy.deepcopy(DefaultPhoSimHeaderMap)
    cat.phoSimHeaderMap['nsnap'] = nsnap
    cat.phoSimHeaderMap['vistime'] = vistime
    if imsim_catalog:
        cat.phoSimHeaderMap['rawSeeing'] = ('rawSeeing', None)
        cat.phoSimHeaderMap['FWHMgeom'] = ('FWHMgeom', None)
        cat.phoSimHeaderMap['FWHMeff'] = ('FWHMeff', None)
    else:
        cat.phoSimHeaderMap['camconfig'] = 1

    with open(outfile, 'w') as output:
        cat.write_header(output)
        if not imsim_catalog:
            output.write('minsource %i\n' % minsource)
            for cat_name in object_catalogs:
                output.write('includeobj %s.gz\n' % cat_name)
    return cat


def get_obs_md(obs_gen, obsHistID, fov=2, dither=True):
    """
    Get the ObservationMetaData object for the specified obsHistID.

    Parameters
    ----------
    obs_gen: lsst.sims.catUtils.utils.ObservationMetaDataGenerator
        Object that reads the opsim db file and generates obs_md objects.
    obsHistID: int
        The ID number of the desired visit.
    fov: float [2]
        Field-of-view angular radius in degrees.  2 degrees will cover
        the LSST focal plane.
    dither: bool [True]
        Flag to apply dithering in the opsim db file.

    Returns
    -------
    lsst.sims.utils.ObservationMetaData object
    """
    obs_md = obs_gen.getObservationMetaData(obsHistID=obsHistID,
                                            boundType='circle',
                                            boundLength=fov)[0]
    if dither:
        obs_md.pointingRA \
            = np.degrees(obs_md.OpsimMetaData['randomDitherFieldPerVisitRA'])
        obs_md.pointingDec \
            = np.degrees(obs_md.OpsimMetaData['randomDitherFieldPerVisitDec'])
        obs_md.OpsimMetaData['rotTelPos'] \
            = obs_md.OpsimMetaData['ditheredRotTelPos']
        obs_md.rotSkyPos \
            = np.degrees(_getRotSkyPos(obs_md._pointingRA, obs_md._pointingDec,
                                       obs_md, obs_md.OpsimMetaData['rotTelPos']))
    return obs_md


def get_instance_catalogs(imsim_catalog=False):
    InstCats = namedtuple('InstCats', ['StarInstCat', 'BrightStarInstCat',
                                       'DESCQACat'])
    if imsim_catalog:
        return InstCats(MaskedPhoSimCatalogPoint_ICRS, BrightStarCatalog_ICRS,
                        PhoSimDESCQA_ICRS)
    return InstCats(MaskedPhoSimCatalogPoint, BrightStarCatalog,
                    PhoSimDESCQA)


class MaskedPhoSimCatalogPoint(PhoSimCatalogPoint):
    disable_proper_motion = False
    min_mag = None
    column_outputs = ['prefix', 'uniqueId', 'raPhoSim', 'decPhoSim',
                      'maskedMagNorm', 'sedFilepath',
                      'redshift', 'gamma1', 'gamma2', 'kappa',
                      'raOffset', 'decOffset',
                      'spatialmodel', 'internalExtinctionModel',
                      'galacticExtinctionModel', 'galacticAv', 'galacticRv']
    protoDc2_half_size = 2.5*np.pi/180.

    @cached
    def get_maskedMagNorm(self):
        raw_norm = self.column_by_name('phoSimMagNorm')
        if self.min_mag is None:
            return raw_norm
        return np.where(raw_norm < self.min_mag, self.min_mag, raw_norm)

    @cached
    def get_inProtoDc2(self):
        ra_values = self.column_by_name('raPhoSim')
        ra = np.where(ra_values < np.pi, ra_values, ra_values - 2.*np.pi)
        dec = self.column_by_name('decPhoSim')
        return np.where((ra > -self.protoDc2_half_size) &
                        (ra < self.protoDc2_half_size) &
                        (dec > -self.protoDc2_half_size) &
                        (dec < self.protoDc2_half_size), 1, None)

    def column_by_name(self, colname):
        if (self.disable_proper_motion and
            colname in ('properMotionRa', 'properMotionDec',
                        'radialVelocity', 'parallax')):
            return np.zeros(len(self.column_by_name('raJ2000')), dtype=np.float)
        return super(MaskedPhoSimCatalogPoint, self).column_by_name(colname)


class BrightStarCatalog(PhoSimCatalogPoint):
    min_mag = None

    @cached
    def get_isBright(self):
        raw_norm = self.column_by_name('phoSimMagNorm')
        return np.where(raw_norm < self.min_mag, raw_norm, None)

class PhoSimDESCQA_ICRS(PhoSimDESCQA):
    catalog_type = 'phoSim_catalog_DESCQA_ICRS'

    column_outputs = ['prefix', 'uniqueId', 'raJ2000', 'decJ2000',
                      'phoSimMagNorm', 'sedFilepath',
                      'redshift', 'gamma1', 'gamma2', 'kappa',
                      'raOffset', 'decOffset',
                      'spatialmodel', 'majorAxis', 'minorAxis',
                      'positionAngle', 'sindex',
                      'internalExtinctionModel', 'internalAv', 'internalRv',
                      'galacticExtinctionModel', 'galacticAv', 'galacticRv',]

    transformations = {'raJ2000': np.degrees,
                       'decJ2000': np.degrees,
                       'positionAngle': np.degrees,
                       'majorAxis': arcsecFromRadians,
                       'minorAxis': arcsecFromRadians}


class MaskedPhoSimCatalogPoint_ICRS(MaskedPhoSimCatalogPoint):
    catalog_type = 'masked_phoSim_catalog_point_ICRS'

    column_outputs = ['prefix', 'uniqueId', 'raJ2000', 'decJ2000',
                      'maskedMagNorm', 'sedFilepath',
                      'redshift', 'gamma1', 'gamma2', 'kappa',
                      'raOffset', 'decOffset',
                      'spatialmodel',
                      'internalExtinctionModel',
                      'galacticExtinctionModel', 'galacticAv', 'galacticRv',]

    transformations = {'raJ2000': np.degrees,
                       'decJ2000': np.degrees}

class BrightStarCatalog_ICRS(BrightStarCatalog):
    catalog_type = 'bright_star_catalog_point_ICRS'

    column_outputs = ['prefix', 'uniqueId', 'raJ2000', 'decJ2000',
                      'phoSimMagNorm', 'sedFilepath',
                      'redshift', 'gamma1', 'gamma2', 'kappa',
                      'raOffset', 'decOffset',
                      'spatialmodel',
                      'internalExtinctionModel',
                      'galacticExtinctionModel', 'galacticAv', 'galacticRv',]

    transformations = {'raJ2000': np.degrees,
                       'decJ2000': np.degrees}
