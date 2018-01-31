import os
import numpy as np
from desc.twinkles.twinklesVariabilityMixins import VariabilityTwinkles
from desc.twinkles import sprinklerCompound
from lsst.sims.catUtils.exampleCatalogDefinitions.phoSimCatalogExamples import PhoSimSpecMap as psmp
from lsst.sims.catalogs.db import CompoundCatalogDBObject
from desc.twinkles import TwinklesCompoundInstanceCatalog, sprinkler
from . import CompoundDESCQAInstanceCatalog, GalaxyCompoundDESCQAObject, PhoSimDESCQA, PhoSimDESCQA_AGN

### This following line is the location you would like to save the SNe SEDs that will be created.
twinkles_sn_sed_dir = 'spectra_files'
twinkles_spec_map = psmp
twinkles_spec_map.subdir_map['(^specFile_)'] = twinkles_sn_sed_dir

_twinkles_defs_file = os.path.join(os.environ['TWINKLES_DIR'], 'data', 'dc2_defs.csv')

_agn_cache_file = os.path.join(os.environ['TWINKLES_DIR'], 'data', 'dc2_agn_cache.csv')
_sne_cache_file = os.path.join(os.environ['TWINKLES_DIR'], 'data', 'dc2_sne_cache.csv')

assert os.path.exists(_agn_cache_file)
assert os.path.exists(_sne_cache_file)

__all__ = ["TwinklesCompoundInstanceCatalog_DC2",
           "sprinklerCompound_DC2", "TwinklesCatalogSersic2D_DC2",
           "TwinklesCatalogZPoint_DC2"]


class sprinklerCompound_DC2(GalaxyCompoundDESCQAObject):
    objid = 'sprinklerCompound_DC2'
    objectTypeId = 166
    cached_sprinkling = True
    agn_cache_file = _agn_cache_file
    sne_cache_file = _sne_cache_file
    defs_file = _twinkles_defs_file

    agn_objid = 'agn_descqa'

    def _final_pass(self, results):

        #Use Sprinkler now
        sp = sprinkler(results, self.mjd, self.specFileMap, density_param=1.0,
                       cached_sprinkling=self.cached_sprinkling,
                       agn_cache_file=self.agn_cache_file,
                       sne_cache_file=self.sne_cache_file,
                       defs_file=self.defs_file)
        results = sp.sprinkle()

        return results

class TwinklesCatalogSersic2D_DC2(PhoSimDESCQA):

    specFileMap = twinkles_spec_map

class TwinklesCatalogZPoint_DC2(PhoSimDESCQA_AGN, VariabilityTwinkles):
    """
    PhoSim Instance Catalog Class for strongly lensed (and therefore time-delayed)
    AGN
    """

    specFileMap = twinkles_spec_map

    catalog_type = 'twinkles_catalog_ZPOINT_DC2'

class TwinklesCompoundInstanceCatalog_DC2(CompoundDESCQAInstanceCatalog):

    use_spec_map = twinkles_spec_map

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

            # These are added to set the field rotation for Proto DC2 catalog
            dbo.field_ra = self._field_ra
            dbo.field_dec = self._field_dec

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

                # These are added here to tell the sprinkler the time so that it
                # can properly create SEDs for supernovae and the spec_map
                # is set to tell the sprinkler where to save those SEDs.
                # There is also a line to set the location of the Proto DC2 AGN DB
                compound_dbo.agn_params_db = self._agn_params_db
                compound_dbo.mjd = self._obs_metadata.mjd.TAI
                compound_dbo.specFileMap = twinkles_spec_map

                self._write_compound(catList, compound_dbo, filename,
                                     chunk_size=chunk_size, write_header=write_header,
                                     write_mode=write_mode)
                write_mode = 'a'
                write_header = False
