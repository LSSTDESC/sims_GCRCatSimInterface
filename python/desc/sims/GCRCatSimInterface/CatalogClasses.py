import os
import numpy as np
from .SedFitter import sed_from_galacticus_mags
from lsst.utils import getPackageDir
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSersic2D
from lsst.sims.catalogs.decorators import cached, compound
from lsst.sims.catUtils.mixins import EBVmixin
from desc.twinkles.twinklesVariabilityMixins import VariabilityTwinkles
from lsst.sims.catUtils.exampleCatalogDefinitions.phoSimCatalogExamples import PhoSimSpecMap as psmp
from lsst.sims.catalogs.definitions import CompoundInstanceCatalog
from lsst.sims.catalogs.db import CompoundCatalogDBObject


twinkles_sn_sed_dir = 'spectra_files'
twinkles_spec_map = psmp
twinkles_spec_map.subdir_map['(^specFile_)'] = twinkles_sn_sed_dir

__all__ = ["PhoSimDESCQA", "TwinklesCompoundInstanceCatalog_DC2"]

#########################################################################
# define a class to write the PhoSim catalog; defining necessary defaults


class TwinklesCompoundInstanceCatalog_DC2(CompoundInstanceCatalog):

    def write_catalog(self, filename, chunk_size=None, write_header=True, write_mode='w'):
        """                                                                                                                                                    
        Write the stored list of InstanceCatalogs to a single ASCII output catalog.                                                                            
        @param [in] filename is the name of the file to be written                                                                                             
        @param [in] chunk_size is an optional parameter telling the CompoundInstanceCatalog                                                                    
        to query the database in manageable chunks (in case returning the whole catalog                                                                        
        takes too much memory)                                                                                                                                 
        @param [in] write_header a boolean specifying whether or not to add a header                                                                           
        to the output catalog (Note: only one header will be written; there will not be                                                                        
        a header for each InstanceCatalog in the CompoundInstanceCatalog; default True)                                                                        
        @param [in] write_mode is 'w' if you want to overwrite the output file or                                                                              
        'a' if you want to append to an existing output file (default: 'w')                                                                                    
        """

        instantiated_ic_list = [None]*len(self._ic_list)

        # first, loop over all of the InstanceCatalog and CatalogDBObject classes, pre-processing                                                              
        # them (i.e. verifying that they have access to all of the columns they need)                                                                   
        for ix, (icClass, dboClass) in enumerate(zip(self._ic_list, self._dbo_list)):
            dbo = dboClass()

            ic = icClass(dbo, obs_metadata=self._obs_metadata)

            # assign all non-private member variables of the CompoundInstanceCatalog                                                 
            # to the instantiated InstanceCatalogs                                                                                                      
            for kk in self.__dict__:
                if kk[0] != '_' and not hasattr(self.__dict__[kk], '__call__'):
                    setattr(ic, kk, self.__dict__[kk])

            for kk in self.__class__.__dict__:
                if kk[0] != '_' and not hasattr(self.__class__.__dict__[kk], '__call__'):
                    setattr(ic, kk, self.__class__.__dict__[kk])

            ic._write_pre_process()
            instantiated_ic_list[ix] = ic

        for row in self._dbObjectGroupList:
            if len(row) == 1:
                ic = instantiated_ic_list[row[0]]
                ic._query_and_write(filename, chunk_size=chunk_size,
                                    write_header=write_header, write_mode=write_mode,
                                    obs_metadata=self._obs_metadata,
                                    constraint=self._constraint)
                write_mode = 'a'
                write_header = False

        default_compound_dbo = None
        if self._compoundDBclass is not None:
            if not hasattr(self._compoundDBclass, '__getitem__'):
                default_compound_dbo = CompoundCatalogDBObject
            else:
                for dbo in self._compoundDBclass:
                    if dbo._table_restriction is None:
                        default_compound_dbo = dbo
                        break

                if default_compound_dbo is None:
                    default_compound_dbo is CompoundCatalogDBObject

        for row in self._dbObjectGroupList:
            if len(row) > 1:
                dbObjClassList = [self._dbo_list[ix] for ix in row]
                catList = [instantiated_ic_list[ix] for ix in row]
                for cat in catList:
                    cat._pre_screen = True

                if self._compoundDBclass is None:
                    compound_dbo = CompoundCatalogDBObject(dbObjClassList)
                elif not hasattr(self._compoundDBclass, '__getitem__'):
                    # if self._compoundDBclass is not a list                                                                                                    
                    try:
                        compound_dbo = self._compoundDBclass(dbObjClassList)
                    except:
                        compound_dbo = default_compound_dbo(dbObjClassList)
                else:
                    compound_dbo = None
                    for candidate in self._compoundDBclass:
                        use_it = True
                        if False in [candidate._table_restriction is not None and
                                     dbo.tableid in candidate._table_restriction
                                     for dbo in dbObjClassList]:

                            use_it = False

                        if use_it:
                            compound_dbo = candidate(dbObjClassList)
                            break

                    if compound_dbo is None:
                        compound_dbo = default_compound_dbo(dbObjClassList)

                compound_dbo.mjd = self._obs_metadata.mjd.TAI
                compound_dbo.specFileMap = twinkles_spec_map

                self._write_compound(catList, compound_dbo, filename,
                                     chunk_size=chunk_size, write_header=write_header,
                                     write_mode=write_mode)
                write_mode = 'a'
                write_header = False

class PhoSimDESCQA(PhoSimCatalogSersic2D, EBVmixin):

    # default values used if the database does not provide information
    default_columns = [('raOffset', 0.0, float), ('decOffset', 0.0, float),
                       ('internalExtinctionModel', 'CCM', str, 3),
                       ('internalAv', 0.1, float),
                       ('internalRv', 3.1, float),
                       ('galacticExtinctionModel', 'CCM', str, 3),
                       ('galacticRv', 3.1, float)]

    # below are defined getter methods used to define CatSim value-added columns
    @cached
    def get_hasDisk(self):
        output = np.where(self.column_by_name('SEDs/diskLuminositiesStellar:SED_9395_583:rest')>0.0, 1.0, None)
        return output

    @cached
    def get_hasBulge(self):
        output = np.where(self.column_by_name('SEDs/spheroidLuminositiesStellar:SED_9395_583:rest')>0.0, 1.0, None)
        return output

    @compound('internalAv', 'internalRv')
    def get_internalDustParams(self):
        if ('hasDisk' in self._cannot_be_null and
            'hasBulge' in self._cannot_be_null):

            raise RuntimeError('\nUnsure whether this is a disk catalog '
                               'or a bulge catalog\n'
                               'self._cannot_be_null %s' % self._cannot_be_null)
        elif 'hasDisk' in self._cannot_be_null:
            lum_type = 'disk'
        elif 'hasBulge' in self._cannot_be_null:
            lum_type = 'spheroid'
        else:
             raise RuntimeError('\nUnsure whether this is a disk catalog '
                               'or a bulge catalog\n'
                               'self._cannot_be_null %s' % self._cannot_be_null)

        b_name = 'otherLuminosities/%sLuminositiesStellar:B:rest' % lum_type
        b_dust_name = 'otherLuminosities/%sLuminositiesStellar:B:rest:dustAtlas' % lum_type

        v_name = 'otherLuminosities/%sLuminositiesStellar:V:rest' % lum_type
        v_dust_name = 'otherLuminosities/%sLuminositiesStellar:V:rest:dustAtlas' % lum_type

        av_list = -2.5*(np.log10(self.column_by_name(v_dust_name)) -
                        np.log10(self.column_by_name(v_name)))

        ebv_list = -2.5*(np.log10(self.column_by_name(b_dust_name)) -
                         np.log10(self.column_by_name(v_dust_name)) -
                         np.log10(self.column_by_name(b_name)) +
                         np.log10(self.column_by_name(v_name)))

        rv_list = av_list/ebv_list

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

    @compound('sedFilename', 'fittedMagNorm')
    def get_fittedSedAndNorm(self):

        if not hasattr(self, '_disk_flux_names'):
            catsim_dir \
                = os.path.join(getPackageDir('sims_GCRCatSimInterface'), 'data')
            catsim_mag_file = os.path.join(catsim_dir, 'CatSimMagGrid.txt')

            if not os.path.exists(catsim_mag_file):
                msg = '\n%s does not exist\n' % catsim_mag_file
                msg += 'Go into the directory %s ' % catsim_dir
                msg += 'and run the script get_sed_mags.py'
                raise RuntimeError(msg)

            with open(catsim_mag_file, 'r') as input_file:
                header = input_file.readlines()[0]
            header = header.strip().split()
            mag_name_list = header[2:-1]
            assert len(mag_name_list) == 30
            bulge_names = []
            disk_names = []
            for name in mag_name_list:
                full_disk_name = 'SEDs/diskLuminositiesStellar:SED_%s:rest' % name
                disk_names.append(full_disk_name)
                full_bulge_name = 'SEDs/spheroidLuminositiesStellar:SED_%s:rest' % name
                bulge_names.append(full_bulge_name)
            self._disk_flux_names = disk_names
            self._bulge_flux_names = bulge_names

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

        sed_names, mag_norms = sed_from_galacticus_mags(mag_array,
                                                        redshift_array)
        return np.array([sed_names, mag_norms])

    def get_phoSimMagNorm(self):
        """
        Need to leave this method here to overload the get_phoSimMagNorm
        in the base PhoSim InstanceCatalog classes
        """
        return self.column_by_name('fittedMagNorm')

