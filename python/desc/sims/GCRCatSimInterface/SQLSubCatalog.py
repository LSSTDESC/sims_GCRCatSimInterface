import numpy as np
import os
import sqlite3
from . import SubCatalogMixin

__all__ = ["SQLSubCatalogMixin"]

class SQLSubCatalogMixin(SubCatalogMixin):

    _table_name = None
    _file_name = None
    _files_written = {}

    def _create_table(self, file_name):
        """
        file_name is the full path to the file where we will
        make the database
        """
        if len(self._current_chunk) is 0:
            return

        dtype_map = {}
        dtype_map[float] = ('float', float)
        dtype_map[np.float] = ('float', float)
        dtype_map[np.float64] = ('float', float)
        dtype_map[np.float32] = ('float', float)
        dtype_map[int] = ('int', int)
        dtype_map[np.int] = ('int', int)
        dtype_map[np.int64] = ('int', int)
        dtype_map[np.int32] = ('int', int)
        dtype_map[np.str_] = ('text', str)
        dtype_map[np.object_] = ('text', str)

        self._type_casting = {}

        with sqlite3.connect(file_name) as conn:
            cursor = conn.cursor()
            creation_cmd = '''CREATE TABLE %s ''' % self._table_name
            creation_cmd += '''('''
            for i_col, name in enumerate(self.iter_column_names()):
                col_type = self.column_by_name(name).dtype.type
                sql_type = dtype_map[col_type][0]
                self._type_casting[name] = dtype_map[col_type][1]
                if i_col>0:
                    creation_cmd += ''', '''
                creation_cmd += '''%s %s''' % (name, sql_type)
            creation_cmd+=''')'''

            cursor.execute(creation_cmd)
            conn.commit()

        if file_name not in self._files_written:
            self._files_written[file_name] = []
        self._files_written[file_name].append(self._table_name)

    def _write_recarray(self, input_recarray, file_handle):
        """
        input_recarray is a recarray of data to be written

        file_handle is the file handle of the main .txt
        InstanceCatalog being written
        """

        if self._table_name is None:
            raise RuntimeError("Cannot call SubCatalogSQLMixin._write_recarray:"
                               "\n_table_name is None")

        if self._file_name is None:
            raise RuntimeError("Cannot call SubCatalogSQLMixin._write_recarray:"
                               "\n_file_name is None")

        self._filter_chunk(input_recarray)

        file_dir = os.path.dirname(file_handle.name)
        full_file_name = os.path.join(file_dir, self._file_name)

        if full_file_name not in self._files_written:
            if os.path.exists(full_file_name):
                os.unlink(full_file_name)

        if (full_file_name not in self._files_written or
            self._table_name not in self._files_written[full_file_name]):

            self._create_table(full_file_name)

        col_dict = {}
        for name in self.iter_column_names():
            arr = self.column_by_name(name)
            if name in self.transformations:
                col_dict[name] = self.transformations[name](arr)
            else:
                col_dict[name] = arr

        if len(self._current_chunk) == 0:
            return

        with sqlite3.connect(full_file_name) as conn:
            insert_cmd = '''INSERT INTO %s ''' % self._table_name
            insert_cmd += '''VALUES('''
            for i_col, name in enumerate(self.iter_column_names()):
                if i_col>0:
                    insert_cmd += ''','''
                insert_cmd += '''?'''
            insert_cmd += ''')'''

            cursor = conn.cursor()
            values = (tuple(self._type_casting[name](col_dict[name][i_obj])
                       for name in self.iter_column_names())
                       for i_obj in range(len(self._current_chunk)))

            cursor.executemany(insert_cmd, values)
            conn.commit()
