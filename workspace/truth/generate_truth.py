from StarTruthModule import write_stars_to_truth
from GalaxyTruthModule import write_galaxies_to_truth

import os
import sqlite3

if __name__ == "__main__":

    db_dir = '/astro/store/pogo3/danielsf/desc_dc2_truth'
    assert os.path.isdir(db_dir)

    param_file = os.path.join(db_dir, 'sprinkled_objects.sqlite')
    assert os.path.isfile(param_file)

    db_file = os.path.join(db_dir, 'proto_dc2_truth_star_gal.db')
    if os.path.isfile(db_file):
        os.unlink(db_file)

    with sqlite3.connect(db_file) as conn:
        cursor = conn.cursor()
        cmd = '''CREATE TABLE truth
              (healpix_2048 int, object_id int, star bool,
               agn bool, sprinkled bool, ra float, dec float,
               redshift float, u float, g float, r float, i float,
               z float, y float)'''

        cursor.execute(cmd)
        conn.commit()

    write_stars_to_truth(output=db_file,
                         n_side=2048,
                         n_procs=20,
                         clobber=False)

    print('wrote stars')

    write_galaxies_to_truth(input=param_file,
                            output=db_file,
                            n_side=2048,
                            n_procs=20,
                            clobber=False)

    print('wrote galaxies')

    with sqlite3.connect(db_file) as conn:
        cursor = conn.cursor()

        cursor.execute('CREATE INDEX obj_id ON truth (object_id)')
        conn.commit()
        print('made object_id index')

        cursor.execute('CREATE INDEX is_star ON truth (is_star)')
        conn.commit()
        print('made is_star index')

        cursor.execute('CREATE INDEX healpix ON truth (healpix_2048)')
        conn.commit()
        print('made healpix index')
