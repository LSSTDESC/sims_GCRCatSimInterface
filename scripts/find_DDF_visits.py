"""
Find visits that intersect DDF region.
"""
import os
import sqlite3
import warnings
from astropy.utils.iers import conf
import numpy as np
import pandas as pd
import lsst.log
import lsst.sphgeom
from lsst.sims.utils import angularSeparation
from lsst.sims.coordUtils import getCornerRaDec, lsst_camera
from desc.sim_utils import DescObsMdGenerator

conf.auto_max_age = None
#local_iers_file = '/home/imSim/data/19-10-30-finals2000A.all'
local_iers_file = '~/dev/imSim/data/19-10-30-finals2000A.all'
conf.iers_auto_url = 'file:' + os.path.abspath(local_iers_file)
conf.auto_download = False

def convex_polygon(corners):
    """
    Return a ConvexPolygon object for the polygon corners represented
    as (longitude, latitude) tuples, in degrees.
    """
    vertices = []
    for corner in corners:
        lonlat = lsst.sphgeom.LonLat.fromDegrees(*corner)
        vertices.append(lsst.sphgeom.UnitVector3d(lonlat))
    return lsst.sphgeom.ConvexPolygon(vertices)

def fp_convex_polygon(detnames, camera, obs_md):
    """
    Return the ConvexPolygons corresponding to outer corners of four
    detector names, assumed to be listed in order of LR, LL, UR, UL
    corners in CCS coordinates.  This ordering is set by the ordering
    of the CCD corners in getCornerRaDec.
    """
    corners = []
    for i, detname in enumerate(detnames):
        corners.append(getCornerRaDec(detname, camera, obs_md)[i])
    return convex_polygon(corners)

class RegionVisitFilter:
    """
    Class to perform selection of visits in an opsim db file based on
    the angular separation between the visit pointing directions and
    the four corners of the region.
    """
    def __init__(self, opsim_db_file, region_corners):
        with sqlite3.connect(opsim_db_file) as conn:
            self.df = pd.read_sql('''select * from summary''', conn)
        self.df['ra'] = np.degrees(self.df['descDitheredRA'].to_numpy())
        self.df['dec'] = np.degrees(self.df['descDitheredDec'].to_numpy())
        for i, corner in enumerate(region_corners):
            self.df[f'sep{i}'] = angularSeparation(self.df['ra'].to_numpy(),
                                                   self.df['dec'].to_numpy(),
                                                   *corner)

    def radius_selection(self, fp_radius=2.0503):
        """
        Query for visits with an offset < fp_radius.
        """
        my_df = self.df.query(f'sep0 < {fp_radius} or sep1 < {fp_radius} or '
                              f'sep2 < {fp_radius} or sep3 < {fp_radius}')
        return my_df

#opsim_db_file = '/home/DC2/minion_1016_desc_dithered_v4_trimmed.db'
opsim_db_file = '/global/cscratch1/sd/jchiang8/desc/Run2.2i/minion_1016_desc_dithered_v4_trimmed.db'
if not os.path.isfile(opsim_db_file):
    raise FileNotFoundError(opsim_db_file)
ddf_corners = ((53.764, -27.533), (52.486, -27.533),
               (52.479, -28.667), (53.771, -28.667))
visit_filter = RegionVisitFilter(opsim_db_file, ddf_corners)

ddf_fiducial = visit_filter.radius_selection(fp_radius=1.76)
ddf_visits = set(ddf_fiducial['obsHistID'])
ddf_inclusive = visit_filter.radius_selection()
candidates = ddf_visits.symmetric_difference(set(ddf_inclusive['obsHistID']))

print(f'{len(ddf_visits)} definite visits')
print(f'processing {len(candidates)} candidate visits:')

lsst.log.setLevel('CameraMapper', lsst.log.WARN)

camera = lsst_camera()
corner_dets_0 = 'R:0,1 S:0,0^R:4,1 S:2,0^R:0,3 S:0,2^R:4,3 S:2,2'.split('^')
corner_dets_1 = 'R:1,0 S:0,0^R:3,0 S:2,0^R:1,4 S:0,2^R:3,4 S:2,2'.split('^')

obs_gen = DescObsMdGenerator(opsim_db_file)

ddf_additional = set()
ddf_polygon = convex_polygon(ddf_corners)
for i, visit in enumerate(candidates):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        obs_md = obs_gen.create(visit)
    for fp_polygon in [fp_convex_polygon(corner_dets_0, camera, obs_md),
                       fp_convex_polygon(corner_dets_1, camera, obs_md)]:
        if ddf_polygon.relate(fp_polygon) != lsst.sphgeom.DISJOINT:
            ddf_additional.add(visit)
            break

# Extract a data frame for the additional visits, and add to the fiducial
# set.
visits = '(' + ','.join([str(_) for _ in ddf_additional]) + ')'
command = f'select * from summary where obsHistID in {visits}'
with sqlite3.connect(opsim_db_file) as conn:
    df = pd.read_sql(command, conn)

ddf_total = pd.concat((ddf_fiducial, df), ignore_index=True, sort=False)

print("total DDF visits:", len(ddf_total))

create_table = """CREATE TABLE IF NOT EXISTS Summary (obsHistID INTEGER, sessionID INTEGER, propID INTEGER, fieldID INTEGER, fieldRA REAL, fieldDec REAL, filter TEXT, expDate INTEGER, expMJD REAL, night INTEGER, visitTime REAL, visitExpTime REAL, finRank REAL, FWHMeff REAL, FWHMgeom REAL, transparency REAL, airmass REAL, vSkyBright REAL, filtSkyBrightness REAL, rotSkyPos REAL, rotTelPos REAL, lst REAL, altitude REAL, azimuth REAL, dist2Moon REAL, solarElong REAL, moonRA REAL, moonDec REAL, moonAlt REAL, moonAZ REAL, moonPhase REAL, sunAlt REAL, sunAz REAL, phaseAngle REAL, rScatter REAL, mieScatter REAL, moonIllum REAL, moonBright REAL, darkBright REAL, rawSeeing REAL, wind REAL, humidity REAL, slewDist REAL, slewTime REAL, fiveSigmaDepth REAL, ditheredRA REAL, ditheredDec REAL, descDitheredDec REAL, descDitheredRA REAL, descDitheredRotTelPos REAL);"""

with sqlite3.connect('static_sims_ddf_visits.db') as conn:
    ddf_total.to_sql('Summary', conn, schema=create_table, index=False)
