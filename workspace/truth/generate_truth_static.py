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
               redshift float,
               u_nodust float, g_nodust float,
               r_nodust float, i_nodust float,
               z_nodust float, y_nodust float,
               u_hostdust float, g_hostdust float,
               r_hostdust float, i_hostdust float,
               z_hostdust float, y_hostdust float,
               u_alldust float, g_alldust float,
               r_alldust float, i_alldust float,
               z_alldust float, y_alldust float)'''

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
                  ('u_nodust', 'observed lsst u magnitude; no dust extinction at all'),
                  ('g_nodust', 'observed lsst g magnitude; no dust extinction at all'),
                  ('r_nodust', 'observed lsst r magnitude; no dust extinction at all'),
                  ('i_nodust', 'observed lsst i magnitude; no dust extinction at all'),
                  ('z_nodust', 'observed lsst z magnitude; no dust extinction at all'),
                  ('y_nodust', 'observed lsst y magnitude; no dust extinction at all'),
                  ('u_hostdust', 'observed lsst u magnitude; only internal dust extinction; no Milky Way dust'),
                  ('g_hostdust', 'observed lsst g magnitude; only internal dust extinction; no Milky Way dust'),
                  ('r_hostdust', 'observed lsst r magnitude; only internal dust extinction; no Milky Way dust'),
                  ('i_hostdust', 'observed lsst i magnitude; only internal dust extinction; no Milky Way dust'),
                  ('z_hostdust', 'observed lsst z magnitude; only internal dust extinction; no Milky Way dust'),
                  ('y_hostdust', 'observed lsst y magnitude; only internal dust extinction; no Milky Way dust'),
                  ('u_alldust', 'observed lsst u magnitude; internal and Milky Way dust extinction applied'),
                  ('g_alldust', 'observed lsst g magnitude; internal and Milky Way dust extinction applied'),
                  ('r_alldust', 'observed lsst r magnitude; internal and Milky Way dust extinction applied'),
                  ('i_alldust', 'observed lsst i magnitude; internal and Milky Way dust extinction applied'),
                  ('z_alldust', 'observed lsst z magnitude; internal and Milky Way dust extinction applied'),
                  ('y_alldust', 'observed lsst y magnitude; internal and Milky Way dust extinction applied')
                  )

        cursor.executemany('INSERT INTO column_descriptions VALUES (?,?)',values)
        conn.commit()

    write_stars_to_truth(output=args.out_file,
                         n_side=2048,
                         n_procs=20,
                         clobber=False)

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
