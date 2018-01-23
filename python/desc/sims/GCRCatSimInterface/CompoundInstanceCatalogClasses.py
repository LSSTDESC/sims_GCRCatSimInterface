from lsst.sims.catalogs.definitions import CompoundInstanceCatalog

class CompoundDESCQAInstanceCatalog(CompoundInstanceCatalog):

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
