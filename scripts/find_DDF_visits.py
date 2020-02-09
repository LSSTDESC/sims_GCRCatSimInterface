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
import lsst.obs.lsst as obs_lsst
import lsst.sphgeom
from lsst.sims.utils import angularSeparation
from lsst.sims.coordUtils import getCornerRaDec
from desc.sim_utils import DescObsMdGenerator

conf.auto_max_age = None
local_iers_file = '/home/imSim/data/19-10-30-finals2000A.all'
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
            self.df = pd.read_sql('''select obsHistID, descDitheredRA,
                                  descDitheredDec from summary''', conn)
        self.df['ra'] = np.degrees(self.df['descDitheredRA'].to_numpy())
        self.df['dec'] = np.degrees(self.df['descDitheredDec'].to_numpy())
        for i, corner in enumerate(region_corners):
            self.df[f'sep{i}'] = angularSeparation(self.df['ra'].to_numpy(),
                                                   self.df['dec'].to_numpy(),
                                                   *corner)

    def get_visits(self, fp_radius=2.0503):
        """
        Query for visits with an offset < fp_radius.
        """
        my_df = self.df.query(f'sep0 < {fp_radius} or sep1 < {fp_radius} or '
                              f'sep2 < {fp_radius} or sep3 < {fp_radius}')
        return set(my_df['obsHistID'])

opsim_db_file = '/home/DC2/minion_1016_desc_dithered_v4_trimmed.db'
ddf_corners = ((53.764, -27.533), (52.486, -27.533),
               (52.479, -28.667), (53.771, -28.667))
visit_filter = RegionVisitFilter(opsim_db_file, ddf_corners)

ddf_visits = visit_filter.get_visits(fp_radius=1.76)
candidates = ddf_visits.symmetric_difference(visit_filter.get_visits())

print(f'{len(ddf_visits)} definite visits')
print(f'processing {len(candidates)} candidate visits:')

lsst.log.setLevel('CameraMapper', lsst.log.WARN)
camera = obs_lsst.LsstCamMapper().camera
obs_gen = DescObsMdGenerator(opsim_db_file)

ddf_polygon = convex_polygon(ddf_corners)
for i, visit in enumerate(candidates):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        obs_md = obs_gen.create(visit)
    for fp_polygon in \
        [fp_convex_polygon('R01_S00 R41_S20 R03_S02 R43_S22'.split(),
                           camera, obs_md),
         fp_convex_polygon('R10_S00 R30_S20 R14_S02 R34_S22'.split(),
                           camera, obs_md)]:
        if ddf_polygon.relate(fp_polygon) != lsst.sphgeom.DISJOINT:
            ddf_visits.add(visit)
            break

print("total DDF visits:", len(ddf_visits))
