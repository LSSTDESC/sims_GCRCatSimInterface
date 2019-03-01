from lsst.sims.catalogs.definitions import CompoundInstanceCatalog

__all__ = ["CompoundDESCQAInstanceCatalog"]

class CompoundDESCQAInstanceCatalog(CompoundInstanceCatalog):

    def __init__(self, instanceCatalogClassList, catalogDBObjectClassList,
                 obs_metadata=None, constraint=None, compoundDBclass = None,
                 field_ra=0., field_dec=0., agn_params_db=None):
        """
        @param [in] instanceCatalogClassList is a list of the InstanceCatalog
        classes to be combined into one output catalog.
        @param [in] catalogDBObjectClassList is a list of the CatalogDBObject
        classes to be associated with the InstanceCatalog classes in
        instanceCatalogClassList.  There should be one CatalogDBObject class
        for each InstanceCatalogClass.
        @param [in] obs_metadata is the ObservationMetaData describing
        the telescope pointing
        @param [in] constraint is an optional SQL constraint to be applied
        to the database query
        @param [in] compoundDBclass is an optional argument specifying what
        CompoundCatalogDBobject class(es) to use to combine InstanceCatalogs
        that query the same table.  This can be either a single
        ComboundCatalogDBObject class, or a list of classes.  The
        CompoundInstanceCatalog will figure out which InstanceCatalog(s) go with
        which CompoundCatalogDBObject class.  If no CompoundCatalogDBObject class
        corresponds to a given group of InstanceCatalogs, then the base
        CompoundCatalogDBObject class will be used.
        Note: compoundDBclass should be a CompoundCatalogDBObject class.
        Not an instantiation of a CompoundCatalogDBObject class.
        """

        self._compoundDBclass = compoundDBclass
        self._obs_metadata = obs_metadata
        self._dbo_list = catalogDBObjectClassList
        self._ic_list = instanceCatalogClassList
        self._constraint = constraint
        self._field_ra = field_ra
        self._field_dec = field_dec
        self._agn_params_db = agn_params_db

        assigned = [False]*len(self._dbo_list)
        self._dbObjectGroupList = []

        for ix in range(len(self._dbo_list)):
            for row in self._dbObjectGroupList:
                if self.areDBObjectsTheSame(self._dbo_list[ix], self._dbo_list[row[0]]):
                    row.append(ix)
                    assigned[ix] = True
                    break

            if not assigned[ix]:
                new_row = [ix]
                for iy in range(ix):
                    if not assigned[iy]:
                        if self.areDBObjectsTheSame(self._dbo_list[ix], self._dbo_list[iy]):
                            new_row.append(iy)

                self._dbObjectGroupList.append(new_row)

    def areDBObjectsTheSame(self, db1, db2):
        """
        Parameters
        ----------
        db1 and db2 are two instantiations of DESCQA CatalogDBObjects.

        Returns
        -------
        True if db1._catalog_id == db2._catalog_id

        False otherwise
        """
        return (db1.yaml_file_name == db2.yaml_file_name and
                db1._cat_cache_suffix == db2._cat_cache_suffix)
