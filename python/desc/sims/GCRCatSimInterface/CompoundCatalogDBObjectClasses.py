from builtins import zip
from builtins import str
from builtins import range
from collections import OrderedDict
from lsst.sims.catalogs.db import _CompoundCatalogDBObject_mixin
from desc.sims.GCRCatSimInterface import DESCQAObject
from desc.sims.GCRCatSimInterface import AGN_postprocessing_mixin

__all__ = ["CompoundDESCQAObject",
           "GalaxyCompoundDESCQAObject"]


class CompoundDESCQAObject(_CompoundCatalogDBObject_mixin, DESCQAObject):
    """
    This is a class for taking several DESCQAObject daughter classes that
    query the same catalog for the same rows (but different
    columns), and combining their queries into one.

    You feed the constructor a list of DESCQAObject daughter classes.  The
    CompoundDESCQAObject verifies that they all do, indeed, query the same
    catalog.  It then constructs its own self.columns member (note
    that CompoundDESCQAObject is a daughter class of DESQAObject) which
    combines all of the requested data.

    When you call query_columns, a recarray will be returned as in a DESCQAObject.
    Note, however, that the names of the columns of the recarray will be modified.
    If the first DESCQAObject class in the list of DESCQAObject classes passed to the
    constructor asks for a column named 'col1', that will be mapped to 'catName_col1'
    where 'catName' is the DESCQAObject's objid member.  'col2' will be mapped to
    'catName_col2', etc.  In cases where the DESCQAObject does not change the name
    of the column, the column will also be returned by its original, un-mangled name.
    """

    def __init__(self, catalogDbObjectClassList, connection=None):
        """
        @param [in] descqaObjectClassList is a list of DESCQAObject
        daughter classes (not instantiations of those classes; the classes
        themselves) that all query the same catalog

        Note: this is a list of classes, not a list of instantiations of those
        classes.  The connection to the catalog is established as soon as
        you instantiate a DESCQAObject daughter class.  To avoid creating
        unnecessary connections, CompoundDESCQAObject will
        read in classes without an active connection and establish its
        own connection in this constructor.  This means that all connection
        parameters must be specified in the class definitions of the classes
        passed into descqaObjectClassList.
        """

        self._dbObjectClassList = catalogDbObjectClassList
        self._nameList = []
        for ix in range(len(self._dbObjectClassList)):
            self._nameList.append(self._dbObjectClassList[ix].objid)

        self._make_columns()
        self.columnMap = OrderedDict([(el[0], el[1] if el[1] else el[0])
                                     for el in self.columns])

        self._validate_input()
        self.objectTypeId = -1  # this is just a placeholder;
                                # the objectTypeId for the classes in
                                # self._dbObjectClassList will actually
                                # be used at runtime

        self.columnMap = dict()
        self._descqaDefaultValues = dict()
        for dbc in self._dbObjectClassList:
            sub_cat_name = dbc.objid
            dbo = dbc()

            for col_name in dbo.descqaDefaultValues:
                prefix_name = '%s_%s'% (sub_cat_name, col_name)
                self._compound_dbo_name_map[prefix_name] = prefix_name
                self._descqaDefaultValues[prefix_name] = dbo.descqaDefaultValues[col_name]

        dbo = self._dbObjectClassList[0]()
        # need to instantiate the first one because sometimes
        # idColKey is not defined until instantiation
        # (see GalaxyTileObj in sims_catUtils/../baseCatalogModels/GalaxyModels.py)

        self.idColKey = dbo.idColKey
        self.yaml_file_name = dbo.yaml_file_name
        self._cat_cache_suffix = dbo._cat_cache_suffix

        super(CompoundDESCQAObject, self).__init__()

    def _make_default_values(self):
        pass

    def _make_column_map(self):
        """
        This dummy method is needed because the call to
        super(CompoundDESCQAObject, self).__init__() in
        this class's __init__ will try to call
        self._make_column_map().  We have already done
        the column map construciton in self.__init__()
        and thus, do not want to do anything more.
        """
        pass

    def _validate_input(self):
        """
        Verify that the DESCQAObjects passed to the constructor
        do, indeed, query the same table of the same catalog.
        """

        dbc0 = self._dbObjectClassList[0]
        for dbc in self._dbObjectClassList:
            if (dbc.yaml_file_name != dbc0.yaml_file_name or
                dbc._cat_cache_suffix != dbc0._cat_cache_suffix):

                raise RuntimeError("Not all DESCQAObject classes "
                                   "passed to CompoundDESCQAObject "
                                   "reference the same catalog")


class GalaxyCompoundDESCQAObject(AGN_postprocessing_mixin, CompoundDESCQAObject):
    agn_params_db = None
    agn_objid = None
