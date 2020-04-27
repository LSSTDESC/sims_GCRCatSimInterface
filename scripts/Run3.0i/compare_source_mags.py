import os
from collections import defaultdict
import sqlite3
import numpy as np
import pandas as pd
import lsst.afw.table as afw_table
import lsst.daf.persistence as dp
import lsst.geom
import desc.sims_ci_pipe as scp


def make_SourceCatalog(df):
    bands = 'ugrizy'
    schema = afw_table.SourceTable.makeMinimalSchema()
    for band in bands:
        schema.addField(f'flux_{band}', type=float, doc=f'{band} flux in nJy')
    src_cat = afw_table.SourceCatalog(schema)
    for iloc in range(len(df)):
        row = df.iloc[iloc]
        new_rec = src_cat.addNew()
        new_rec.set('id', int(row['id']))
        new_rec.set('coord_ra', lsst.geom.Angle(row.ra, lsst.geom.degrees))
        new_rec.set('coord_dec', lsst.geom.Angle(row.dec, lsst.geom.degrees))
        for band in bands:
            colname = f'flux_{band}'
            new_rec.set(colname, row[colname])
    return src_cat


def get_point_sources(butler, visit, flux_type='base_PsfFlux'):
    datarefs = butler.subset('src', visit=visit)
    dfs = []
    for dataref in list(datarefs):
        print(dataref.dataId)
        src = scp.get_point_sources(dataref.get('src'))
        calexp = dataref.get('calexp')
        calib = calexp.getPhotoCalib()
        flux, fluxerr = calib.instFluxToNanojansky(src, flux_type).transpose()

        dfs.append(pd.DataFrame(data=dict(ra=np.degrees(src.get('coord_ra')),
                                          dec=np.degrees(src.get('coord_dec')),
                                          flux=flux, fluxerr=fluxerr)))
    return pd.concat(dfs)


def match_meas_fluxes(butler, visit, star_truth_summary_file,
                      flux_type='base_PsfFlux', max_offset=0.1):
    flux_col = f'{flux_type}_instFlux'
    conn = sqlite3.connect(star_truth_summary_file)
    radius = lsst.geom.Angle(max_offset, lsst.geom.arcseconds)
    dfs = []
    datarefs = butler.subset('src', visit=visit)
    for i, dataref in enumerate(list(datarefs)):
        print(i)
        calib = dataref.get('calexp').getPhotoCalib()
        src = scp.get_point_sources(dataref.get('src'))
        ras = np.degrees(src.get('coord_ra'))
        decs = np.degrees(src.get('coord_dec'))
        ra_min, ra_max = min(ras), max(ras)
        dec_min, dec_max = min(decs), max(decs)
        query = f'''select * from truth_summary where
                    {ra_min} <= ra and ra <= {ra_max} and
                    {dec_min} <= dec and dec <= {dec_max}'''
        truth_df = pd.read_sql(query, conn)
        truth_cat = make_SourceCatalog(truth_df)

        matches = afw_table.matchRaDec(truth_cat, src, radius)
        num_matches = len(matches)

        ids = np.zeros(num_matches, dtype=np.int)
        offsets = np.zeros(num_matches, dtype=np.float)
        true_fluxes = np.zeros(num_matches, dtype=np.float)
        meas_fluxes = np.zeros(num_matches, dtype=np.float)
        meas_fluxerrs = np.zeros(num_matches, dtype=np.float)

        for i, match in enumerate(matches):
            ids[i] = match.first['id']
            offsets[i] = np.degrees(match.distance)*3600*1000.
            true_fluxes[i] = match.first[f'flux_{band}']
            meas_fluxes[i] = calib.instFluxToNanojansky(match.second[flux_col])
            meas_fluxerrs[i] \
                = calib.instFluxToNanojansky(match.second[flux_col + 'Err'])

        dfs.append(pd.DataFrame(data=dict(id=ids, offset=offsets,
                                          true_flux=true_fluxes,
                                          meas_flux=meas_fluxes,
                                          meas_fluxerr=meas_fluxerrs)))
    df = pd.concat(dfs)
    return df


def compute_delta_fluxes(stars_db_file, ids, mjd):
    import desc.sims_truthcatalog as stc
    lc_factory = stc.StellarLightCurveFactory(stars_db_file)
    mjds = [mjd]
    delta_fluxes = defaultdict(dict)
    for obj_id in ids:
        dm, m0 = lc_factory.create(obj_id, mjds)
        for band in dm:
            delta_fluxes[obj_id][band] \
                = (10.**((8.9 - (m0[band] + dm[band]))/2.5)
                   - 10.**((8.9 - m0[band])/2.5))*1e9
    return delta_fluxes


if __name__ == '__main__':
    star_truth_summary_file = '/global/cscratch1/sd/descim/star_truth/star_truth_summary.db'

    repo = 'repo_lensed_sne'
    butler = dp.Butler(repo)

    visit = 709692
    band = 'i'

    outfile = f'src_truth_match_v{visit}-{band}.pkl'
    if not os.path.isfile(outfile):
        df = match_meas_fluxes(butler, visit, star_truth_summary_file)
        df.to_pickle(outfile)
    else:
        df = pd.read_pickle(outfile)

