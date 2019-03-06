import sqlite3
import numpy as np

from lsst.sims.utils import ObservationMetaData
from lsst.sims.utils import htmModule as htm
from lsst.sims.catUtils.baseCatalogModels import StarObj

in_db = StarObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                port=1433, driver='mssql+pymssql')

# full DC2 field of view
dc2_obs = ObservationMetaData(pointingRA=55.064, pointingDec=-29.783,
                              boundType='circle', boundLength=23.0)

# file to write
out_db_name = 'dc2_stellar_db.db'

with sqlite3.connect(out_db_name) as out_conn:
    out_c = out_conn.cursor()
    creation_query = '''CREATE TABLE stars (
                     simobjid int, htmid_6 int, ra real, decl real,
                     gal_l real, gal_b real, magNorm real,
                     mura real, mudecl real, parallax real,
                     ebv real, radialVelocity real, varParamStr text,
                     sedFilename text,
                     umag real, gmag real, rmag real, imag real,
                     zmag real, ymag real)'''
    out_c.execute(creation_query)
    out_conn.commit()

    colnames = ['simobjid', 'ra', 'decl', 'gal_l', 'gal_b',
                'magNorm', 'mura', 'mudecl', 'parallax', 'ebv',
                'radialVelocity', 'varParamStr', 'sedFilename',
                'umag', 'gmag', 'rmag', 'imag', 'zmag', 'ymag']

    data_iter = in_db.query_columns(colnames=colnames, obs_metadata=dc2_obs,
                                    chunk_size=100000)

    for data_chunk in data_iter:

        # spatial indexing to help with queries
        htmid = htm.findHtmid(data_chunk['ra'], data_chunk['decl'], 6)

        vals = ((int(d['simobjid']), int(hh),
                 d['ra'], d['decl'], d['gal_l'], d['gal_b'],
                 d['magNorm'], d['mura'], d['mudecl'], d['parallax'],
                 d['ebv'], d['radialVelocity'], d['varParamStr'],
                 d['sedFilename'], d['umag'], d['gmag'],
                 d['rmag'], d['imag'], d['zmag'], d['ymag'])
                for d, hh in zip(data_chunk, htmid))
        out_c.executemany('''INSERT INTO stars VALUES
                          (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', vals)
        out_conn.commit()

    out_c.execute('''CREATE INDEX htmid_6_idx ON stars (htmid_6)''')
    out_conn.commit()
