import numpy as np
import os
import sqlite3

__all__ = ["SubCatalogSQLMixin"]

class SubCatalogSQLMixin(object):

    _table_name = None
    _file_name = None
    _files_written = set()
    _tables_created = set()

    def _create_table(self, file_name, local_recarray):
        """
        file_name is the full path to the file where we will
        make the database

        local_recarray contains the data to be written
        """
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

        self._type_casting = {}

        with sqlite3.connect(file_name) as conn:
            cursor = conn.cursor()
            creation_cmd = '''CREATE TABLE %s ''' % self._table_name
            creation_cmd += '''('''
            for i_col, name in enumerate(local_recarray.dtype.names):
                col_type = local_recarray.dtype[name]
                sql_type = dtype_map[col_type][0]
                self._type_casting[name] = dtype_map[col_type][1]
                if i_col>0:
                    creation_cmd += ''', '''
                creation_cmd += '''%s %s''' % (name, sql_type)
            creation_cmd+=''')'''

            cursor.execute()
            conn.commit()

        self._files_written.add(file_name)
        self._tables_created.add(self._table_name)

    def _write_recarray(self, local_recarray, file_handle):
        """
        local_recarray is a recarray of data to be written

        file_handle is the file handle of the main .txt
        InstanceCatalog being written
        """

        if self._table_name is None:
            raise RuntimeError("Cannot call SubCatalogSQLMixin._write_recarray:"
                               "\n_table_name is None")

        if self._file_name is None:
            raise RuntimeError("Cannot call SubCatalogSQLMixin._write_recarray:"
                               "\n_file_name is None")

        file_dir = os.path.dirname(file_handle.name)
        full_file_name = os.path.join(file_dir, self._file_name)

        if full_file_name not in self._files_written:
            if os.path.exists(full_file_name):
                os.unlink(full_file_name)

        if self._table_name not in self._tables_created:
            self._create_table(full_file_name, local_recarray)

        with sqlite3.connect(full_file_name) as conn:
            insert_cmd = '''INSERT INTO %s ''' % self._table_name
            insert_cmd += '''VALUES('''
            for i_col in range(len(local_recarray.dtype.names)):
                if i_col>0:
                    insert_cmd += ''','''
                insert_cmd += '''?'''
            insert_cmd += ''')'''

            cursor = conn.cursor()
            values = ((self._type_casting[name](local_recarray[name][i_obj])
                       for name in local_recarray.dtype.names)
                       for i_obj in range(len(local_recarray)))
            cursor.executemany(insert_cmd, values)
            conn.commit()
