"""
This is the script that actually generates the
truth catalog of static sources based on the database
of parameters generated in generate_truth_params.

It took about 1.5 hours to run for protoDC2.
"""

from desc.sims.GCRCatSimInterface.StarTruthModule import write_stars_to_truth
from desc.sims.GCRCatSimInterface.GalaxyTruthModule import write_galaxies_to_truth

import os
import sqlite3
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--param_file', type=str,
                        default=os.path.join('/astro/store/pogo3',
                                             'danielsf/desc_dc2_truth',
                                             'sprinkled_objects.sqlite'))
    parser.add_argument('--out_file', type=str,
                       default=os.path.join('/astro/store/pogo3',
                                            'danielsf/desc_dc2_truth',
                                            'proto_dc2_truth_star_gal.db'))

    args = parser.parse_args()

    assert os.path.isfile(args.param_file)

    if os.path.isfile(args.out_file):
        os.unlink(args.out_file)

    with sqlite3.connect(args.out_file) as conn:
        cursor = conn.cursor()
        cmd = '''CREATE TABLE truth
              (healpix_2048 int, object_id int, star int,
               agn int, sprinkled int, ra float, dec float,
               redshift float, u float, g float, r float, i float,
               z float, y float, du_mw float, dg_mw float,
               dr_mw float, di_mw float, dz_mw float,
               du_internal float, dg_internal float,
               dr_internal float, di_internal float,
               dz_internal_float, dy_internal_float)'''

        cursor.execute(cmd)
        conn.commit()

        cmd = '''CREATE TABLE column_descriptions
              (name text, description text)'''
        cursor.execute(cmd)

        values = (('healpix_2048', 'healpixel containing the object (nside=2048; nested)'),
                  ('object_id', 'an int uniquely identifying objects (can collide between stars, galaxies, and sprinkled objects)'),
                  ('star', 'an int; ==1 if a star; ==0 if not'),
                  ('agn', 'an int; ==1 if galaxy has an AGN; ==0 if not'),
                  ('sprinkled', 'an int; ==1 if object added by the sprinkler; ==0 if not'),
                  ('ra', 'in degrees'),
                  ('dec', 'in degrees'),
                  ('redshift', 'cosmological only'),
                  ('u', 'observed lsst u magnitude'),
                  ('g', 'observed lsst g magnitude'),
                  ('r', 'observed lsst r magnitude'),
                  ('i', 'observed lsst i magnitude'),
                  ('z', 'observed lsst z magnitude'),
                  ('y', 'observed lsst y magnitude'),
                  ('du_mw', 'magnitudes of extinction in lsst u due to Milky Way dust'),
                  ('dg_mw', 'magnitudes of extinction in lsst g due to Milky Way dust'),
                  ('dr_mw', 'magnitudes of extinction in lsst r due to Milky Way dust'),
                  ('di_mw', 'magnitudes of extinction in lsst i due to Milky Way dust'),
                  ('dz_mw', 'magnitudes of extinction in lsst z due to Milky Way dust'),
                  ('dy_mw', 'magnitudes of extinction in lsst y due to Milky Way dust'),
                  ('du_internal', 'magnitudes of extinction in lsst u due to internal dust'),
                  ('dg_internal', 'magnitudes of extinction in lsst g due to internal dust'),
                  ('dr_internal', 'magnitudes of extinction in lsst r due to internal dust'),
                  ('di_internal', 'magnitudes of extinction in lsst i due to internal dust'),
                  ('dz_internal', 'magnitudes of extinction in lsst z due to internal dust'),
                  ('dy_internal', 'magnitudes of extinction in lsst y due to internal dust'))


        cursor.executemany('INSERT INTO column_descriptions VALUES (?,?)',values)
        conn.commit()

    #write_stars_to_truth(output=args.out_file,
    #                     n_side=2048,
    #                     n_procs=20,
    #                     clobber=False)

    print('wrote stars')

    write_galaxies_to_truth(input_db=args.param_file,
                            output=args.out_file,
                            n_side=2048,
                            n_procs=20,
                            clobber=False)

    print('wrote galaxies')

    with sqlite3.connect(args.out_file) as conn:
        cursor = conn.cursor()

        cursor.execute('CREATE INDEX obj_id ON truth (object_id)')
        conn.commit()
        print('made object_id index')

        cursor.execute('CREATE INDEX is_star ON truth (star)')
        conn.commit()
        print('made is_star index')

        cursor.execute('CREATE INDEX is_sprinkled ON truth (sprinkled)')
        conn.commit()
        print('made is_sprinkled index')

        cursor.execute('CREATE INDEX healpix ON truth (healpix_2048)')
        conn.commit()
        print('made healpix index')
