import os
import sqlite3
import gzip
import numpy as np

import lsst.sims.utils.htmModule as htm
from lsst.sims.utils import angularSeparation
import lsst.sims.photUtils as photUtils
from lsst.sims.catUtils.supernovae import SNObject

import time
import argparse


def validate_sne(cat_dir, obsid, fov_deg=2.1,
                 scatter_file=None, vanishing_file=None):
    """
    Parameters
    ----------
    cat_dir is the parent dir of $obsid

    obsid is the obsHistID of the pointing (an int)

    fov_deg is the radius of the field of view in degrees
    """

    project_dir = os.path.join('/global/projecta/projectdirs',
                               'lsst/groups/SSim/DC2/cosmoDC2_v1.1.4')

    sne_db_name = os.path.join(project_dir, 'sne_cosmoDC2_v1.1.4_MS_DDF.db')
    if not os.path.isfile(sne_db_name):
        raise RuntimeError("\n%s\nis not a file" % sne_db_name)

    instcat_dir = os.path.join(cat_dir, '%.8d' % obsid)
    if not os.path.isdir(instcat_dir):
        raise RuntimeError("\n%s\nis not a dir" % instcat_dir)

    sne_name = os.path.join(instcat_dir, 'sne_cat_%d.txt.gz' % obsid)
    if not os.path.isfile(sne_name):
        raise RuntimeError("\n%s\nis not a file" % sne_name)

    phosim_name = os.path.join(instcat_dir, 'phosim_cat_%d.txt' % obsid)
    if not os.path.isfile(phosim_name):
        raise RuntimeError("\n%s\nis not a file" % phosim_name)

    sed_dir = os.path.join(instcat_dir, "Dynamic")
    if not os.path.isdir(sed_dir):
        raise RuntimeError("\n%s\nis not a dir" % sed_dir)

    opsim_db = os.path.join("/global/projecta/projectdirs",
                            "lsst/groups/SSim/DC2",
                            "minion_1016_desc_dithered_v4_sfd.db")
    if not os.path.isfile(opsim_db):
        raise RuntimeError("\n%s\nis not a file" % opsim_db)

    band_from_int = 'ugrizy'
    bandpass = None
    with open(phosim_name, 'r') as in_file:
        for line in in_file:
            params = line.strip().split()
            if params[0] == 'filter':
                bandpass = band_from_int[int(params[1])]
                break
    if bandpass is None:
        raise RuntimeError("Failed to read in bandpass")

    bp_dict = photUtils.BandpassDict.loadTotalBandpassesFromFiles()

    with sqlite3.connect(opsim_db) as conn:
        c = conn.cursor()
        r = c.execute("SELECT descDitheredRA, descDitheredDec, expMJD FROM Summary "
                      "WHERE obsHistID==%d" % obsid).fetchall()
        pointing_ra = np.degrees(r[0][0])
        pointing_dec = np.degrees(r[0][1])
        expmjd = float(r[0][2])

    with sqlite3.connect(sne_db_name) as conn:
        c = conn.cursor()
        query = "SELECT snid_in, snra_in, sndec_in, "
        query += "c_in, mB, t0_in, x0_in, x1_in, z_in "
        query += "FROM sne_params WHERE ("

        where_clause = None
        half_space = htm.halfSpaceFromRaDec(pointing_ra, pointing_dec, fov_deg)
        bounds = half_space.findAllTrixels(6)

        for bound in bounds:
            if where_clause is not None:
                where_clause += " OR ("
            else:
                where_clause = "("

            if bound[0] == bound[1]:
                where_clause += "htmid_level_6==%d" % bound[0]
            else:
                where_clause += "htmid_level_6>=%d AND htmid_level_6<=%d" % (bound[0],bound[1])

            where_clause += ")"

        query += where_clause+")"

        sn_params = c.execute(query).fetchall()
    sn_ra = np.array([sn[1] for sn in sn_params])
    sn_dec = np.array([sn[2] for sn in sn_params])

    sn_d = angularSeparation(pointing_ra, pointing_dec,
                             sn_ra, sn_dec)

    valid = np.where(sn_d<=fov_deg)

    sn_ra_dec = {}
    sn_param_dict = {}
    for i_sn in valid[0]:
        sn = sn_params[i_sn]
        sn_param_dict[sn[0]] = sn[3:]
        sn_ra_dec[sn[0]] = sn[1:3]

    sne_in_instcat = set()
    d_mag_max = -1.0
    with gzip.open(sne_name, 'rb') as sne_cat_file:
        for sne_line in sne_cat_file:
            instcat_params = sne_line.strip().split(b' ')
            sne_id = instcat_params[1].decode('utf-8')
            if sne_id not in sn_param_dict:
                try:
                    sne_id_int = int(sne_id)
                    assert sne_id_int<1.5e10
                except (ValueError, AssertionError):
                    raise RuntimeError("\n%s\nnot in SNe db" % sne_id)

            sne_in_instcat.add(sne_id)

            control_params = sn_param_dict[sne_id]
            sed_name = os.path.join(instcat_dir, instcat_params[5].decode('utf-8'))
            if not os.path.isfile(sed_name):
                raise RuntimeError("\n%s\nis not a file" % sed_name)

            sed = photUtils.Sed()
            sed.readSED_flambda(sed_name, cache_sed=False)

            fnorm = photUtils.getImsimFluxNorm(sed, float(instcat_params[4]))
            sed.multiplyFluxNorm(fnorm)

            sed.redshiftSED(float(instcat_params[6]), dimming=True)
            a_x, b_x = sed.setupCCM_ab()
            A_v = float(instcat_params[15])
            R_v = float(instcat_params[16])
            EBV = A_v/R_v
            sed.addDust(a_x, b_x, A_v=A_v, R_v=R_v)

            sn_object = SNObject()
            sn_object.set_MWebv(EBV)
            sn_object.set(c=control_params[0], x1=control_params[4],
                          x0=control_params[3], t0=control_params[2],
                          z=control_params[5])

            sn_mag = sn_object.catsimBandMag(bp_dict[bandpass], expmjd)
            instcat_flux = sed.calcFlux(bp_dict[bandpass])
            instcat_mag = sed.magFromFlux(instcat_flux)
            d_mag = (sn_mag-instcat_mag)
            if not np.isfinite(d_mag):
                if np.isfinite(sn_mag) or np.isfinite(instcat_mag):
                    msg = '%s\ngave magnitudes %e %e (diff %e)' % (sne_id,
                                                                   sn_mag,
                                                                   instcat_mag,
                                                                   d_mag)
                    print(msg)

            if scatter_file is not None:
                scatter_file.write("%e %e %e\n" % (sn_mag, instcat_mag, d_mag))

            d_mag = np.abs(d_mag)
            if d_mag>d_mag_max:
                d_mag_max = d_mag
                #print('dmag %e -- %e (%e)' %
                #      (d_mag, instcat_mag,
                #       control_params[5]-float(instcat_params[6])))

    msg = '\n%s\n' % sne_name
    ct_missing = 0
    for sn_id in sn_param_dict:
        if sn_id in sne_in_instcat:
            continue
        control_params = sn_param_dict[sn_id]

        sn_object = SNObject()
        sn_object.set_MWebv(0.0)
        cc = control_params[0]
        x1 = control_params[4]
        x0 = control_params[3]
        t0 = control_params[2]
        zz = control_params[5]
        sn_object.set(c=cc, x1=x1, x0=x0, t0=t0, z=zz)

        sn_mag = sn_object.catsimBandMag(bp_dict[bandpass], expmjd)
        if vanishing_file is not None and np.isfinite(sn_mag):
            vanishing_file.write('%d %s %s %e ' % (obsid, bandpass, sn_id, sn_mag))
            vanishing_file.write('%e %e %e %.4f %e %.4f %e\n' %
                                 (cc,x0,x1,t0,zz,expmjd,expmjd-t0))

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--cat_dir', type=str, default=None,
                        help='parent dir of $obs')
    parser.add_argument('--obs', type=int, default=None,
                        help='obsHistID of the pointing')
    parser.add_argument('--fov_deg', type=float, default=2.1,
                        help='radius of field of view in degrees')

    args = parser.parse_args()

    t_start = time.time()
    validate_sne(args.cat_dir, args.obs, fov_deg=args.fov_deg)
    print('that took %e seconds' % (time.time()-t_start))
