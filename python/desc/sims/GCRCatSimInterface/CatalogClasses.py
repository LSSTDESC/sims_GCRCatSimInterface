import os
import re
import numpy as np
import healpy
import h5py
import copy
from lsst.utils import getPackageDir
from desc.sims.GCRCatSimInterface import _DESCQAObject_metadata
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.decorators import cached
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSersic2D
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogZPoint
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSN
from lsst.sims.utils import angularSeparation
from lsst.sims.catUtils.mixins import VariabilityAGN
from lsst.sims.catalogs.decorators import cached, compound
from lsst.sims.catUtils.mixins import EBVmixin


__all__ = ["PhoSimDESCQA", "PhoSimDESCQA_AGN", "DC2PhosimCatalogSN",
           "SubCatalogMixin", "SprinklerTruthCatMixin", "TruthPhoSimDESCQA",
           "TruthPhoSimDESCQA_AGN"]

#########################################################################
# define a class to write the PhoSim catalog; defining necessary defaults


class SubCatalogMixin(object):
    """
    This mixin provides a way to write parallel catalogs from
    a CompoundInstanceCatalog.  It supplants the _write_recarray
    method, which CompundInstanceCatalog calls, and replaces it
    with something that will write a separate truth catalog.
    """

    _subcat_file_handle = None

    subcat_prefix = None  # prefix prepended to main InstanceCatalog file name
    subcat_suffix = None  # suffix appended to main InstanceCatalog file name

    # This boolean will keep track of whether or not this
    # truth catalog has been written to yet.  If it has,
    # it will be opened in mode 'a'; if not, it will be
    # opened in mode 'w'
    _subcat_cat_written = False

    # so that we don't have to rework the CompoundInstanceCatalog
    # API
    _write_subcat_header = False

    # The list below *will* be shared among instantiations
    # as a safeguard against accidentally opening the
    # same SubCatalog in write mode twice
    _list_of_opened_subcats = set()

    @cached
    def get_sprinkling_switch(self):
        is_sprinkled = self.column_by_name('is_sprinkled')
        return np.where(is_sprinkled==1, 1, None)

    def _get_subcat_name(self, file_handle):
        """
        file_handle points to the main InstanceCatalog that
        the CompoundInstanceCatalog is trying to write

        Returns a string containing the name of the subcatalog
        that this subcatalog instantiation is writing.
        """

        file_dir = os.path.dirname(file_handle.name)
        instcat_name = os.path.basename(file_handle.name)
        subcat_file_name = instcat_name

        if self.subcat_prefix is None and self.subcat_suffix is None:
            raise RuntimeError("Trying to write SubCatalog without either "
                               "a subcat_prefix or a subcat_suffix. This "
                               "could cause you to overwrite existing files")

        if self.subcat_prefix is not None:
            subcat_file_name = self.subcat_prefix + subcat_file_name
        if self.subcat_suffix is not None:
            subcat_file_name += self.subcat_suffix

        subcat_name = os.path.join(file_dir,
                                   subcat_file_name)

        assert subcat_name != file_handle.name
        return subcat_name

    def _write_recarray(self, local_recarray, file_handle):
        """
        local_recarray is a recarray of the data to be written

        file_handle points to the main InstanceCatalog that
        the CompoundInstanceCatalog is trying to write
        """
        if self._subcat_file_handle is None:
            subcat_name = self._get_subcat_name(file_handle)
            if not self._subcat_cat_written:
                write_mode = 'w'
                if subcat_name in self._list_of_opened_subcats:
                    raise RuntimeError("Trying to create SubCatalog\n"
                                       + "%s\n" % subcat_name
                                       + "which was already created")
            else:
                write_mode = 'a'
            self._subcat_file_handle = open(subcat_name, write_mode)
            self._subcat_cat_written = True
            self._list_of_opened_subcats.add(subcat_name)

            if write_mode == 'w' and self._write_subcat_header:
                # call InstanceCatalog.write_header to avoid calling
                # the PhoSim catalog write_header (which will require
                # a phoSimHeaderMap)
                InstanceCatalog.write_header(self, self._subcat_file_handle)

        InstanceCatalog._write_recarray(self, local_recarray,
                                        self._subcat_file_handle)

        self._subcat_file_handle.flush()


class SprinklerTruthCatMixin(SubCatalogMixin):
    """
    A sub-class of the SubCatalogMixin specifically for generating truth
    catalogs for sprinkled objects
    """
    column_outputs = ['uniqueId', 'galaxy_id', 'raJ2000', 'decJ2000',
                      'sedFilepath', 'phoSimMagNorm',
                      'redshift', 'isPoint']

    cannot_be_null = ['sprinkling_switch']

    _write_subcat_header = True


class DC2PhosimCatalogSN(PhoSimCatalogSN):
    """
    Modification of the PhoSimCatalogSN mixin to provide shorter sedFileNames
    by leaving out the parts of the directory name. Also fix name changes from
    gamma to shear.
    """
    def get_reasonableMagNorm(self):
        with np.errstate(invalid='ignore', divide='ignore'):
            mn = self.column_by_name('phoSimMagNorm')
            return np.where(mn<500.0, mn, np.NaN)

    def get_uniqueId(self):
        return self.column_by_name(self.refIdCol)

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

    cannot_be_null = ['x0', 't0', 'z', 'reasonableMagNorm']

    default_columns = [('gamma1', 0., float), ('gamma2', 0., float), ('kappa', 0., float),
                       ('raOffset', 0., float), ('decOffset', 0., float),
                       ('galacticAv', 0.1, float), ('galacticRv', 3.1, float),
                       ('galacticExtinctionModel', 'CCM', (str, 3)),
                       ('internalExtinctionModel', 'none', (str, 4)), ('internalAv', 0., float),
                       ('internalRv', 3.1, float), ('shear1', 0., float), ('shear2', 0., float)]


class PhoSimDESCQA(PhoSimCatalogSersic2D, EBVmixin):

    # directory where the SED lookup tables reside
    sed_lookup_dir = None

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

        super(PhoSimDESCQA, self).__init__(*args, **kwargs)

    def get_component_type(self):
        """
        returns 'disk' if this is a catalog disks;
        returns 'bulge' if this is a catalog of bulges
        """
        if not hasattr(self, '_lumtype'):

            if ('hasDisk' in self._cannot_be_null and
                'hasBulge' in self._cannot_be_null):

                raise RuntimeError('\nUnsure whether this is a disk catalog '
                                   'or a bulge catalog\n'
                                   'self._cannot_be_null %s' % self._cannot_be_null)
            elif 'hasDisk' in self._cannot_be_null:
                self._lum_type = 'disk'
            elif 'hasKnots' in self._cannot_be_null:
                self._lum_type = 'knots'
            elif 'hasBulge' in self._cannot_be_null:
                self._lum_type = 'bulge'
            else:
                raise RuntimeError('\nUnsure whether this is a disk catalog '
                                   'or a bulge catalog\n'
                                   'self._cannot_be_null %s' % self._cannot_be_null)

        return self._lum_type

    # below are defined getter methods used to define CatSim value-added columns
    @cached
    def get_hasDisk(self):
        output = np.where(self.column_by_name('stellar_mass_disk')>0.0,
                          1.0, np.NaN)
        return output

    @cached
    def get_hasKnots(self):
        return self.column_by_name('hasDisk')

    @cached
    def get_hasBulge(self):
        output = np.where(self.column_by_name('stellar_mass_bulge')>0.0, 1.0, np.NaN)
        return output

    def _cache_sed_lookup(self, healpix_list, component_type, bandpass):
        """
        Load the SED lookup table information for the healpixels specified
        in healpix_list.

        component_type is either 'disk' or 'bulge'

        bandpass is one of 'ugrizy'

        Return as a dict pointing to numpy arrays with
        keys:

        galaxy_id
        *_sed
        *_magnorm
        *_av
        *_rv

        where * stands for either 'disk' or 'bulge'
        """

        if component_type != 'disk' and component_type != 'bulge':
            raise RuntimeError("Do not know what component this is: %s" % component_type)


        assert os.path.isdir(self.sed_lookup_dir)

        file_root = 'sed_fit'
        bp_to_int = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'y':5}

        out_dict = {}

        raw_out_dict = {}
        raw_out_dict['galaxy_id'] = []
        raw_out_dict['sed_idx'] = []
        raw_out_dict['magnorm'] = []
        raw_out_dict['av'] = []
        raw_out_dict['rv'] = []

        for hp in healpix_list:
            file_name = os.path.join(self.sed_lookup_dir, '%s_%d.h5' % (file_root, hp))
            with h5py.File(file_name, 'r') as data:
                if not hasattr(self, '_sed_lookup_names'):
                    self._sed_lookup_names = np.copy(data['sed_names']).astype(str)
                    self._sed_lookup_names_bytes = np.copy(data['sed_names'])
                else:
                    np.testing.assert_array_equal(data['sed_names'],
                                                  self._sed_lookup_names_bytes)

                dd = angularSeparation(self.obs_metadata.pointingRA,
                                       self.obs_metadata.pointingDec,
                                       data['ra'].value, data['dec'].value)
                to_keep = np.where(dd<self.obs_metadata.boundLength+0.1)

                raw_out_dict['galaxy_id'].append(data['galaxy_id'].value[to_keep])
                raw_out_dict['sed_idx'].append(data['%s_sed' % component_type].value[to_keep])
                raw_out_dict['magnorm'].append(data['%s_magnorm' % component_type].value[bp_to_int[bandpass]][to_keep])
                raw_out_dict['av'].append(data['%s_av' % component_type].value[to_keep])
                raw_out_dict['rv'].append(data['%s_rv' % component_type].value[to_keep])



        out_dict['galaxy_id'] = np.concatenate(raw_out_dict.pop('galaxy_id'))

        out_dict['%s_sed_idx'
                 % component_type] = np.concatenate(raw_out_dict.pop('sed_idx'))

        out_dict['%s_%s_magnorm'
                 % (component_type,
                    bandpass)] = np.concatenate(raw_out_dict.pop('magnorm'))

        out_dict['%s_av'
                 % component_type] = np.concatenate(raw_out_dict.pop('av'))

        out_dict['%s_rv'
                 % component_type] = np.concatenate(raw_out_dict.pop('rv'))

        # so that we can use numpy search sorted
        sorted_dex = np.argsort(out_dict['galaxy_id'])
        for k in out_dict.keys():
            out_dict[k] = out_dict[k][sorted_dex]

        return out_dict

    @compound('sedFilename_idx', 'magNorm_fitted',
              'internalAv_fitted', 'internalRv_fitted')
    def get_fittedSedAndNorm(self):

        _knots_cutoff_i_mag = 27.0

        component_type = self.get_component_type()

        self.column_by_name('raJ2000')
        self.column_by_name('decJ2000')
        galaxy_id = self.column_by_name('galaxy_id')

        if not hasattr(self, '_knots_available'):
            self._knots_available = False
            if 'knots_flux_ratio' in self.db_obj._catalog.list_all_quantities(include_native=True):
                self._knots_available = True

        if self._knots_available:

            if not hasattr(self, '_sprinkled_gid'):
                twinkles_dir = os.path.join(os.environ['TWINKLES_DIR'], 'data')
                agn_name = os.path.join(twinkles_dir, 'cosmoDC2_v1.1.4_agn_cache.csv')
                sne_name = os.path.join(twinkles_dir, 'cosmoDC2_v1.1.4_sne_cache.csv')
                sprinkled_gid = []
                for file_name in (agn_name, sne_name):
                    with open(file_name, 'r') as in_file:
                        for line in in_file:
                            if line.startswith('galtileid'):
                                continue
                            params = line.strip().split(',')
                            sprinkled_gid.append(int(params[0]))
                self._sprinkled_gid = np.array(sprinkled_gid)

            lsst_i_mag = self.column_by_name('mag_true_i_lsst')
            knots_ratio = self.column_by_name('knots_flux_ratio')
            is_sprinkled = np.in1d(galaxy_id, self._sprinkled_gid,
                                   assume_unique=True)

            knots_ratio = np.where(~is_sprinkled,
                                   knots_ratio, 0.0)

        if component_type == 'knots' and not self._knots_available:
            raise RuntimeError("You are trying to simulate knots "
                               "but there are no knots in your "
                               "extragalactic catalog")

        if len(galaxy_id) == 0:
            return np.array([[],[], [], []])

        if hasattr(self.db_obj, '_loaded_healpixel'):
            healpix_list = np.array([self.db_obj._loaded_healpixel])
        elif (hasattr(self, 'filter_on_healpix') and
              self.filter_on_healpix is True and
              'loaded_healpixel' in _DESCQAObject_metadata):

            healpix_list = np.array([_DESCQAObject_metadata['loaded_healpixel']])
        else:
            ra_rad = self.obs_metadata._pointingRA
            dec_rad = self.obs_metadata._pointingDec
            vv = np.array([np.cos(dec_rad)*np.cos(ra_rad),
                           np.cos(dec_rad)*np.sin(ra_rad),
                           np.sin(dec_rad)])
            radius_rad = self.obs_metadata._boundLength
            healpix_list = np.sort(healpy.query_disc(32, vv, radius_rad,
                                                     inclusive=True, nest=False))

        if component_type == 'knots':
            cache_component_type = 'disk'
        else:
            cache_component_type = component_type

        if (not hasattr(self, '_sed_lookup_cache') or
            not np.array_equal(healpix_list, self._sed_lookup_healpix) or
            not component_type == self._sed_lookup_component_type or
            not self.obs_metadata.bandpass == self._sed_lookup_bandpass):

            # discard any existing cache
            if hasattr(self, '_sed_lookup_cache'):
                del self._sed_lookup_cache

            self._sed_lookup_cache = self._cache_sed_lookup(healpix_list,
                                                            cache_component_type,
                                                            self.obs_metadata.bandpass)
            self._sed_lookup_healpix = np.copy(healpix_list)
            self._sed_lookup_bandpass = self.obs_metadata.bandpass
            self._sed_lookup_component_type = component_type

        idx = np.searchsorted(self._sed_lookup_cache['galaxy_id'],
                              galaxy_id)

        # The sprinkler will add some galaxy_id that do not map to the SED
        # lookup cache (the sprinkler will add SEDs for these galaxies, so
        # we do not need to worry about fitting an SED to them).  These
        # galaxies can be identified because their galaxy_id values will
        # be larger than any galaxy in the extragalactic catalog, thus,
        # the idx values assigned by np.searchsorted will == len(_sed_lookup_cache).
        # We now remove those galaxies from the fitting.
        valid_gal = np.where(idx<len(self._sed_lookup_cache['galaxy_id']))
        idx = idx[valid_gal]

        np.testing.assert_array_equal(self._sed_lookup_cache['galaxy_id'][idx],
                                      galaxy_id[valid_gal])

        n_gal = len(galaxy_id)
        sed_idx = -1*np.ones(n_gal, dtype=int)
        mag_norms = np.NaN*np.ones(n_gal, dtype=float)
        av = np.NaN*np.ones(n_gal, dtype=float)
        rv = np.NaN*np.ones(n_gal, dtype=float)

        sed_idx[valid_gal] = self._sed_lookup_cache['%s_sed_idx' % cache_component_type][idx]
        mag_norms[valid_gal] = self._sed_lookup_cache['%s_%s_magnorm' % (cache_component_type, self.obs_metadata.bandpass)][idx]
        av[valid_gal] = self._sed_lookup_cache['%s_av' % cache_component_type][idx]
        rv[valid_gal] = self._sed_lookup_cache['%s_rv' % cache_component_type][idx]

        with np.errstate(invalid='ignore', divide='ignore'):
            if component_type != 'bulge' and self._knots_available:
                if component_type == 'disk':
                    d_mag = np.where(lsst_i_mag<=_knots_cutoff_i_mag,
                                     -2.5*np.log10(1.0-knots_ratio), 0.0)
                elif component_type == 'knots':
                    d_mag = np.where(lsst_i_mag<=_knots_cutoff_i_mag,
                                     -2.5*np.log10(knots_ratio), np.NaN)
                else:
                    raise RuntimeError("Not sure how to handle d_mag for component %s" % component_type)

                mag_norms += d_mag

        return np.array([sed_idx, mag_norms, av, rv], dtype=object)

    @cached
    def get_magNorm(self):
        if self.get_component_type() == 'bulge':
            magnorm_name = 'bulgeMagNorm'
        else:
            magnorm_name = 'diskMagNorm'
        with np.errstate(invalid='ignore', divide='ignore'):
            raw_magnorm = self.column_by_name(magnorm_name)
            fitted_magnorm = self.column_by_name('magNorm_fitted')
            preliminary_output=np.where(np.isnan(raw_magnorm), fitted_magnorm, raw_magnorm)
            preliminary_output = np.array(preliminary_output).astype(float)
            return np.where(preliminary_output<998.0, preliminary_output, np.NaN)

    @cached
    def get_sedFilepath(self):
        if self.get_component_type() == 'bulge':
            sed_name = 'bulgeSedFilename'
        else:
            sed_name = 'diskSedFilename'
        raw_filename = self.column_by_name(sed_name)
        sed_idx = self.column_by_name('sedFilename_idx')
        if len(sed_idx)==0:
            return np.array([])

        fitted_filename = np.array([self._sed_lookup_names[ii]
                                    if ii>=0 else 'None'
                                    for ii in sed_idx])

        return np.where(np.char.find(raw_filename.astype('str'), 'None')==0,
                        fitted_filename, raw_filename)

    @cached
    def get_internalRv(self):
        if self.get_component_type() == 'bulge':
            rv_name = 'bulgeInternalRv'
        else:
            rv_name = 'diskInternalRv'
        raw_rv = self.column_by_name(rv_name)
        fitted_rv = self.column_by_name('internalRv_fitted')
        return np.where(np.isnan(raw_rv), fitted_rv, raw_rv)


    @cached
    def get_internalAv(self):
        if self.get_component_type() == 'bulge':
            av_name = 'bulgeInternalAv'
        else:
            av_name = 'diskInternalAv'
        raw_av = self.column_by_name(av_name)
        fitted_av = self.column_by_name('internalAv_fitted')
        return np.where(np.isnan(raw_av), fitted_av, raw_av)

    def get_phoSimMagNorm(self):
        """
        Need to leave this method here to overload the get_phoSimMagNorm
        in the base PhoSim InstanceCatalog classes
        """
        self.column_by_name('is_sprinkled')
        return self.column_by_name('magNorm')


class TruthPhoSimDESCQA(SprinklerTruthCatMixin, PhoSimDESCQA):
    pass

class PhoSimDESCQA_AGN(PhoSimCatalogZPoint, EBVmixin, VariabilityAGN):

    column_outputs = ['prefix', 'uniqueId', 'raPhoSim', 'decPhoSim', 'magNormFiltered', 'sedFilepath',
                      'redshift', 'gamma1', 'gamma2', 'kappa', 'raOffset', 'decOffset',
                      'spatialmodel', 'internalExtinctionModel',
                      'galacticExtinctionModel', 'galacticAv', 'galacticRv']


    cannot_be_null = ['magNormFiltered']

    @cached
    def get_magNorm(self):
        gid = self.column_by_name('galaxy_id')
        bogey = 2563036006
        dx = np.where(gid==bogey)
        for bp in 'ugrizy':
            dmag = self.column_by_name('delta_lsst_%s' % bp)
            if len(dx[0])>0:
                print('delta_%s: %e' % (bp,dmag[dx]))
        if len(dx[0])>0:
            print('mjd %e' % (self.obs_metadata.mjd.TAI))
        return self.column_by_name('agnMagNorm')

    @cached
    def get_sedFilename(self):
        return self.column_by_name('agnSedFilename')

    @cached
    def get_magNormFiltered(self):
        mm = self.column_by_name('phoSimMagNorm')
        return np.clip(mm, 10.0, None)

    @cached
    def get_prefix(self):
        self.column_by_name('is_sprinkled')
        chunkiter = range(len(self.column_by_name(self.refIdCol)))
        return np.array(['object' for i in chunkiter], dtype=(str, 6))


class TruthPhoSimDESCQA_AGN(SprinklerTruthCatMixin, PhoSimDESCQA_AGN):
    pass
