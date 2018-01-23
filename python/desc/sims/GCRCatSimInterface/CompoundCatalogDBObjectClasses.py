from lsst.sims.catalogs.db import CompoundCatalogDBObject

__all__ = ["CompoundDESCQACatalogDBObject"]

class CompoundDESCQACatalogDBObject(CompoundCatalogDBObject):

    def _make_dbTypeMap(self):
        pass

    def _make_dbDefaultValues(self):
        pass

    def _validate_input(self):
        """
        Verify that CatalogDBObject classes passed to the constructor
        do, indeed, all query the same catalog.
        """

        dbc0 = self._dbObjectClassList[0]

        acceptable = True

        for dbc in self._dbObjectClassList:
            if dbc.yaml_file_name != dbc0.yaml_file_name:
                acceptable = False
            if dbc._cat_cache_suffix != dbc0._cat_cache_suffix:
                acceptable = False

        if not acceptable:
            raise RuntimeError("The CatalogDBObject classes passed to "
                               "CompoundDESCQACatalogDBObject do not all "
                               "query the same catalog")
