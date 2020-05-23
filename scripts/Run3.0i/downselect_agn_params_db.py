import sqlite3
import pandas as pd

agn_params_file = '/global/cscratch1/sd/jchiang8/desc/sims_GCRCatSimInterface/work/2020-02-18/agn_cosmoDC2_v1.1.4.db'
schema = '''CREATE TABLE agn_params
                       (galaxy_id int, htmid_8 int, magNorm real,
                        redshift real, M_i real, ra real, dec real,
                        varParamStr text);
CREATE INDEX htmid ON agn_params (htmid_8);
'''

ra_min, ra_max = 52.479, 53.771
dec_min, dec_max = -28.667, -27.533

query = f'''select * from agn_params where {ra_min} <= ra and ra <= {ra_max}
            and {dec_min} <= dec and dec <= {dec_max}'''

with sqlite3.connect(agn_params_file) as conn:
    df = pd.read_sql(query, conn)

outfile = 'agn_cosmoDC2_v1.1.4_ddf.db'
with sqlite3.connect(outfile) as conn:
    df.to_sql('agn_params', conn, schema=schema, index=False)
