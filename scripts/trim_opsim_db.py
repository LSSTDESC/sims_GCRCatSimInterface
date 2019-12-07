import os
import sqlite3

wfd_visit_list = os.path.join(os.environ['SIMS_GCRCATSIMINTERFACE_DIR'],
                              'workspace', 'run2.1', 'data',
                              'master_obshistid_list.txt')
with open(wfd_visit_list) as fd:
    wfd_visits = [line.strip() for line in fd]

# SQL command to the summary table.
create_table = """CREATE TABLE IF NOT EXISTS Summary (obsHistID INTEGER, sessionID INTEGER, propID INTEGER, fieldID INTEGER, fieldRA REAL, fieldDec REAL, filter TEXT, expDate INTEGER, expMJD REAL, night INTEGER, visitTime REAL, visitExpTime REAL, finRank REAL, FWHMeff REAL, FWHMgeom REAL, transparency REAL, airmass REAL, vSkyBright REAL, filtSkyBrightness REAL, rotSkyPos REAL, rotTelPos REAL, lst REAL, altitude REAL, azimuth REAL, dist2Moon REAL, solarElong REAL, moonRA REAL, moonDec REAL, moonAlt REAL, moonAZ REAL, moonPhase REAL, sunAlt REAL, sunAz REAL, phaseAngle REAL, rScatter REAL, mieScatter REAL, moonIllum REAL, moonBright REAL, darkBright REAL, rawSeeing REAL, wind REAL, humidity REAL, slewDist REAL, slewTime REAL, fiveSigmaDepth REAL, ditheredRA REAL, ditheredDec REAL, descDitheredDec REAL, descDitheredRA REAL, descDitheredRotTelPos REAL);"""

src_conn = sqlite3.connect('/home/DC2/minion_1016_desc_dithered_v4.db')
outfile = 'minion_1016_desc_dithered_v4_trimmed.db'
if os.path.isfile(outfile):
    os.remove(outfile)
dest_conn = sqlite3.connect(outfile)

dest_curs = dest_conn.cursor()
dest_curs.execute(create_table)
dest_conn.commit()

visits = "(" + ", ".join(wfd_visits) + ")"
get_rows = f"""select * from summary where obshistid in {visits}"""
src_curs = src_conn.cursor()
src_curs.execute(get_rows)
rows = src_curs.fetchall()

slots = ','.join(['?' for i in range(len(rows[0]))])
insert_rows = f'insert into summary values ({slots})'
dest_curs.executemany(insert_rows, rows)
dest_conn.commit()

create_indexes = """CREATE INDEX fieldID_idx ON Summary(fieldID);
CREATE INDEX expMJD_idx ON Summary(expMJD);
CREATE INDEX filter_idx ON Summary(filter);
CREATE INDEX fieldRA_idx ON Summary(fieldRA);
CREATE INDEX fieldDec_idx ON Summary(fieldDec);
CREATE INDEX fieldRADec_idx ON Summary(fieldRA, fieldDec);
CREATE INDEX night_idx ON Summary(night);
CREATE INDEX propID_idx ON Summary(propID);
CREATE INDEX ditheredRA_idx ON Summary(ditheredRA);
CREATE INDEX ditheredDec_idx ON Summary(ditheredDec);
CREATE INDEX ditheredRADec_idx ON Summary(ditheredRA, ditheredDec);
CREATE INDEX filter_propID_idx ON Summary(filter, propID);
""".split('\n')
for command in create_indexes:
    dest_curs.execute(command)
dest_conn.commit()
