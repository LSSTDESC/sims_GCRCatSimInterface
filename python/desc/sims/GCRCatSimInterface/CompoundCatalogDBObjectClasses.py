from builtins import zip
from builtins import str
from builtins import range
from desc.sims.GCRCatSimInterface import DESCQAObject

__all__ = ["CompoundDESCQAObject"]


class CompoundDESCQAObject(DESCQAObject):
    """
    This is a class for taking several CatalogDBObject daughter classes that
    query the same table of the same database for the same rows (but different
    columns; note that the columns can be transformed by the CatalogDBObjects'
    self.columns member), and combining their queries into one.

    You feed the constructor a list of CatalogDBObject daughter classes.  The
    CompoundCatalogDBObject verifies that they all do, indeed, query the same table
    of the same database.  It then constructs its own self.columns member (note
    that CompoundCatalogDBObject is a daughter class of CatalogDBObject) which
    combines all of the requested data.

    When you call query_columns, a recarray will be returned as in a CatalogDBObject.
    Note, however, that the names of the columns of the recarray will be modified.
    If the first CatalogDBObject in the list of CatalogDBObjects passed to the constructor
    asks for a column named 'col1', that will be mapped to 'catName_col1' where 'catName'
    is the CatalogDBObject's objid member.  'col2' will be mapped to 'catName_col2', etc.
    In cases where the CatalogDBObject does not change the name of the column, the column
    will also be returned by its original, un-mangled name.

    In cases where a custom query_columns method must be implemented, this class
    can be sub-classed and the custom method added as a member method.  In that
    case, the _table_restriction member variable should be set to a list of table
    names corresponding to the tables for which this class was designed.  An
    exception will be raised if the user tries to use the CompoundCatalogDBObject
    class to query tables for which it was not written.  _table_restriction defaults
    to None, which means that the class is for use with any table.
    """

    def __init__(self, catalogDbObjectClassList, connection=None):
        """
        @param [in] catalogDbObjectClassList is a list of CatalogDBObject
        daughter classes (not instantiations of those classes; the classes
        themselves) that all query the same database table

        Note: this is a list of classes, not a list of instantiations of those
        classes.  The connection to the database is established as soon as
        you instantiate a CatalogDBObject daughter class.  To avoid creating
        unnecessary database connections, CompoundCatalogDBObject will
        read in classes without an active connection and establish its
        own connection in this constructor.  This means that all connection
        parameters must be specified in the class definitions of the classes
        passed into catalogDbObjectClassList.

        @param [in] connection is an optional instantiation of DBConnection
        representing an active connection to the database required by
        this CompoundCatalogDBObject (prevents the CompoundCatalogDBObject
        from opening a redundant connection)
        """

        self._dbObjectClassList = catalogDbObjectClassList
        self._validate_input()
        self.objectTypeId = 119

        self.columnMap = dict()
        for dbc in self._dbObjectClassList:
            sub_cat_name = dbc.objid
            dbo = dbc()
            for col_name in dbo.columnMap:
                self.columnMap[sub_cat_name+'_'+col_name] = dbo.columnMap[col_name]

        dbo = self._dbObjectClassList[0]()
        # need to instantiate the first one because sometimes
        # idColKey is not defined until instantiation
        # (see GalaxyTileObj in sims_catUtils/../baseCatalogModels/GalaxyModels.py)

        self.idColKey = dbo.idColKey
        self.yaml_file_name = dbo.yaml_file_name
        self._cat_cache_suffix = dbo._cat_cache_suffix

        super(CompoundDESCQAObject, self).__init__()

    def _make_column_map(self):
        pass

    def _validate_input(self):
        """
        Verify that the CatalogDBObjects passed to the constructor
        do, indeed, query the same table of the same database.
        """

        dbc0 = self._dbObjectClassList[0]
        for dbc in self._dbObjectClassList:
            if (dbc.yaml_file_name != dbc0.yaml_file_name or
                dbc._cat_cache_suffix != dbc0._cat_cache_suffix):

                raise RuntimeError("Not all DESCQAObject classes "
                                   "passed to CompoundDESCQAObject "
                                   "reference the same catalog")
