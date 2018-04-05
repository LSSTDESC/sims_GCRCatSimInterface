import os
import re
import numpy as np
import copy
from .SedFitter import sed_from_galacticus_mags
from .SedFitter import sed_filter_names_from_catalog
from lsst.utils import getPackageDir
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSersic2D
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogZPoint
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSN
from lsst.sims.catUtils.mixins import VariabilityAGN
from lsst.sims.catalogs.decorators import cached, compound
from lsst.sims.catUtils.mixins import EBVmixin


__all__ = ["PhoSimDESCQA", "PhoSimDESCQA_AGN", "DC2PhosimCatalogSN"]

#########################################################################
# define a class to write the PhoSim catalog; defining necessary defaults


class DC2PhosimCatalogSN(PhoSimCatalogSN):
    """
    Modification of the PhoSimCatalogSN mixin to provide shorter sedFileNames
    by leaving out the parts of the directory name. Also fix name changes from
    gamma to shear.
    """
    def get_sedFilepath(self):
        return self.column_by_name('TsedFilepath')

    def get_shorterFileNames(self):
        """
        Method to truncate filenames for transient
        spectra written out by phosim.

        .. note: the variable sep needs to be in
        `self.sn_sedfile_prefix` before writing out
        a phosim catalog.
        """
        fnames = self.column_by_name('sedFilepath')
        sep = 'Dynamic/specFileSN_'
        split_names = []
        for fname in fnames:
            if 'None' not in fname:
                fname = sep + fname.split(sep)[-1]
            else:
                fname = 'None'
            split_names.append(fname)
        return np.array(split_names)

    # column_outputs = PhoSimCatalogSN.column_outputs
    # column_outputs[PhoSimCatalogSN.column_outputs.index('sedFilepath')] = \
    #    'shorterFileNames'
    column_outputs = ['prefix', 'uniqueId', 'raPhoSim', 'decPhoSim',
                      'phoSimMagNorm', 'shorterFileNames', 'redshift',
                      'shear1', 'shear2', 'kappa', 'raOffset', 'decOffset',
                      'spatialmodel', 'internalExtinctionModel',
                      'galacticExtinctionModel', 'galacticAv', 'galacticRv']

    cannot_be_null = ['x0', 't0', 'z', 'shorterFileNames']

    default_columns = [('gamma1', 0., float), ('gamma2', 0., float), ('kappa', 0., float),
                       ('raOffset', 0., float), ('decOffset', 0., float),
                       ('galacticAv', 0.1, float), ('galacticRv', 3.1, float),
                       ('galacticExtinctionModel', 'CCM', (str, 3)),
                       ('internalExtinctionModel', 'none', (str, 4)), ('internalAv', 0., float),
                       ('internalRv', 3.1, float), ('shear1', 0., float), ('shear2', 0., float)]


class PhoSimDESCQA(PhoSimCatalogSersic2D, EBVmixin):

    # default values used if the database does not provide information
    default_columns = [('raOffset', 0.0, float), ('decOffset', 0.0, float),
                       ('internalExtinctionModel', 'CCM', str, 3),
                       ('internalAv', 0.1, float),
                       ('internalRv', 3.1, float),
                       ('galacticExtinctionModel', 'CCM', str, 3),
                       ('galacticRv', 3.1, float)]

    cannot_be_null = ['magNorm']

    def __init__(self, *args, **kwargs):
        # Update the spatial model if knots are requested, for knots, the sersic
        # parameter actually contains the number of knots
        if 'cannot_be_null' in kwargs.keys():
            if 'hasKnots' in kwargs['cannot_be_null']:
                self.catalog_type = 'phoSim_catalog_KNOTS'
                self.spatialModel = 'knots'
                if 'hasDisk' not in kwargs['cannot_be_null']:
                    kwargs['cannot_be_null'].append('hasDisk')

        super(PhoSimDESCQA, self).__init__(*args, **kwargs)

    # below are defined getter methods used to define CatSim value-added columns
    @cached
    def get_hasDisk(self):
        output = np.where(self.column_by_name('stellar_mass_disk')>0.0, 1.0, None)
        return output

    @cached
    def get_hasKnots(self):
        return self.column_by_name('hasDisk')

    @cached
    def get_hasBulge(self):
        output = np.where(self.column_by_name('stellar_mass_bulge')>0.0, 1.0, None)
        return output

    @compound('internalAv_fitted', 'internalRv_fitted')
    def get_internalDustParams(self):
        if ('hasDisk' in self._cannot_be_null and
            'hasBulge' in self._cannot_be_null):

            raise RuntimeError('\nUnsure whether this is a disk catalog '
                               'or a bulge catalog\n'
                               'self._cannot_be_null %s' % self._cannot_be_null)
        elif 'hasDisk' in self._cannot_be_null:
            lum_type = 'disk'
        elif 'hasBulge' in self._cannot_be_null:
            lum_type = 'bulge'
        else:
             raise RuntimeError('\nUnsure whether this is a disk catalog '
                               'or a bulge catalog\n'
                               'self._cannot_be_null %s' % self._cannot_be_null)

        av_name = 'A_v_%s' % lum_type
        if av_name not in self._all_available_columns:
            av_name = 'A_v'
        av_list = copy.copy(self.column_by_name(av_name))

        rv_name = 'R_v_%s' % lum_type
        if rv_name not in self._all_available_columns:
            rv_name = 'R_v'
        rv_list = copy.copy(self.column_by_name(rv_name))

        # this is a hack to replace anomalous values of dust extinction
        # with more reasonable values
        if not hasattr(self, '_dust_rng'):
            self._dust_rng = np.random.RandomState(182314)

        offensive_av = np.where(np.logical_or(av_list<0.001, av_list>3.1))

        av_list[offensive_av] = self._dust_rng.random_sample(len(offensive_av[0]))*3.1+0.001

        offensive_rv = np.where(np.logical_or(np.isnan(rv_list),
                                np.logical_or(rv_list<1.0, rv_list>5.0)))
        rv_list[offensive_rv] = self._dust_rng.random_sample(len(offensive_rv[0]))*4.0+1.0

        return np.array([av_list, rv_list])


    @compound('sedFilename_fitted', 'magNorm_fitted')
    def get_fittedSedAndNorm(self):

        if not hasattr(self, '_disk_flux_names'):

            f_params = sed_filter_names_from_catalog(self.db_obj._catalog)

            np.testing.assert_array_almost_equal(f_params['disk']['wav_min'],
                                                 f_params['bulge']['wav_min'],
                                                 decimal=10)

            np.testing.assert_array_almost_equal(f_params['disk']['wav_width'],
                                                 f_params['bulge']['wav_width'],
                                                 decimal=10)

            self._disk_flux_names = f_params['disk']['filter_name']
            self._bulge_flux_names = f_params['bulge']['filter_name']
            self._sed_wav_min = f_params['disk']['wav_min']
            self._sed_wav_width = f_params['disk']['wav_width']

        if 'hasBulge' in self._cannot_be_null and 'hasDisk' in self._cannot_be_null:
            raise RuntimeError('\nUnsure whether this is a disk catalog or a bulge catalog.\n'
                               'Both appear to be in self._cannot_be_null.\n'
                               'self._cannot_be_null: %s' % self._cannot_be_null)
        elif 'hasBulge' in self._cannot_be_null:
            flux_names = self._bulge_flux_names
        elif 'hasDisk' in self._cannot_be_null:
            flux_names = self._disk_flux_names
        else:
            raise RuntimeError('\nUnsure whether this is a disk catalog or a bluge catalog.\n'
                               'Neither appear to be in self._cannot_be_null.\n'
                               'self._cannot_be_null: %s' % self._cannot_be_null)

        mag_array = np.array([-2.5*np.log10(self.column_by_name(name))
                              for name in flux_names])

        redshift_array = self.column_by_name('true_redshift')

        if len(redshift_array) == 0:
            return np.array([[], []])

        H0 = self.db_obj._catalog.cosmology.H0.value
        Om0 = self.db_obj._catalog.cosmology.Om0

        (sed_names,
         mag_norms) = sed_from_galacticus_mags(mag_array,
                                               redshift_array,
                                               H0, Om0,
                                               wav_min=self._sed_wav_min,
                                               wav_width=self._sed_wav_width)

        return np.array([sed_names, mag_norms])

    @cached
    def get_magNorm(self):
        raw_magnorm = self.column_by_name('magNorm_dc2')
        fitted_magnorm = self.column_by_name('magNorm_fitted')
        preliminary_output=np.where(np.isnan(raw_magnorm), fitted_magnorm, raw_magnorm)
        preliminary_output = np.array(preliminary_output).astype(float)
        return np.where(preliminary_output<998.0, preliminary_output, np.NaN)

    @cached
    def get_sedFilename(self):
        raw_filename = self.column_by_name('sedFilename_dc2')
        fitted_filename = self.column_by_name('sedFilename_fitted')
        return np.where(np.char.find(raw_filename.astype('str'), 'None')==0,
                        fitted_filename, raw_filename)

    @cached
    def get_internalRv(self):
        raw_rv = self.column_by_name('internalRv_dc2')
        fitted_rv = self.column_by_name('internalRv_fitted')
        return np.where(np.isnan(raw_rv), fitted_rv, raw_rv)


    @cached
    def get_internalAv(self):
        raw_av = self.column_by_name('internalAv_dc2')
        fitted_av = self.column_by_name('internalAv_fitted')
        return np.where(np.isnan(raw_av), fitted_av, raw_av)

    def get_phoSimMagNorm(self):
        """
        Need to leave this method here to overload the get_phoSimMagNorm
        in the base PhoSim InstanceCatalog classes
        """
        return self.column_by_name('magNorm')

class PhoSimDESCQA_AGN(PhoSimCatalogZPoint, EBVmixin, VariabilityAGN):

    cannot_be_null = ['sedFilepath', 'magNorm']
