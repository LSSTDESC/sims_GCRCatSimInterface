"""
Requires snsims, healpy, pandas and lsst_sims

Example command line calls for mDDF and Main Survey:
python write_healpixel.py --data_root /Users/rbiswas/data/DESC/DC2 --survey MDDF
python write_healpixel.py --data_root /Users/rbiswas/data/DESC/DC2 --survey MS --file_root gals_564_ra_dec.hdf --healpixelId 564

Assumes that the mddf field is entirely in a single healpixel for mDDF SN
For running mDDF, assumes that `sn_564_MS.csv` has been written out to disk in pwd
"""
from __future__ import absolute_import, division, print_function
import numpy as np
import pandas as pd
import healpy as hp
from argparse import ArgumentParser
from joinSNCat import DC2SN
from joinSNCat import DC2 as dc2Cosmo
import snsims
import os

import GCRCatalogs

if __name__ == '__main__':
    

    parser = ArgumentParser()
    parser.add_argument('--out_dir', help='directory to store cvs files in',
                        default='./')
    parser.add_argument('--data_root', help='path to directory where the galaxy files are located, defaults="./"',
                        default='./')
    parser.add_argument('--file_root', help='file name of the galaxy hdf5 file without the healpixel, eg. "gals_{0}_ra_dec.hdf" to which this will default to',
                        default='gals_{0}_ra_dec.hdf')
    parser.add_argument('--healpixelId', type=int, help='healpixel ID used', default=None)
    parser.add_argument('--randomSeedOffset', type=int, help='integer to offset random seeds by, None offsets by healpixelID', default=None)
    parser.add_argument('--survey', type=str, help='main survey or DDF, MS|mDDF', default='MS')
    parser.add_argument('--zmax', type=float, help='max redshift of SN, defaults to None, if None, zmax=1.0 for MS, 1.4 for DDF', default=None)

    args = parser.parse_args()

    print(snsims.__version__)
    print(args)

    assert os.path.isdir(args.out_dir)
    NSIDE = 32
    zmax = args.zmax

    if args.healpixelId is not None:
        healpix_list = [args.healpixelId]
    elif args.survey == 'mDDF':
        ra_rad = np.radians(53.125)
        dec_rad = np.radians(-28.10)
        vv = np.array([np.cos(dec_rad)*np.cos(ra_rad),
                       np.cos(dec_rad)*np.sin(ra_rad),
                       np.sin(dec_rad)])
        healpix_list = hp.query_disc(32, vv, np.radians(0.6),
                                     nest=False, inclusive=True)
    else:
        cf = GCRCatalogs.get_catalog_config('cosmoDC2_v1.1.4_image')
        healpix_list = cf['healpix_pixels']

    for healpixelId in healpix_list:
        fname = os.path.join(args.data_root, args.file_root)
        fname = fname.format(healpixelId)
        survey = args.survey

        if args.randomSeedOffset is None:
            randomSeedOffset = healpixelId
        else:
            randomSeedOffset = args.randomSeedOffset

        veto_ms_gals = False

        if survey.lower() == 'mddf':
            healpixelSN_fname = os.path.join(os.environ['SCRATCH'],
                                             'cosmoDC2_v1.1.4_sne',
                                             'sne_csv',
                                             'sn_%d_MS.csv' % healpixelId)
            assert os.path.isfile(healpixelSN_fname)
            print('reading veto data from %s' % healpixelSN_fname)
            #healpixelSN_fname = 'sn_564_MS.csv'
            randomSeedOffset = 20000+healpixelId
            veto_ms_gals = True
            # This is the size of region in square degrees
            area = 1.28 #  dsinthetadtheta = np.cos(np.radians(90 + 27.53)) -  np.cos(np.radians(90 + 28.667))
            # dphi = np.radians(53.76 - 52.48)
            if zmax is None:
                zmax = 1.4
        else:
            # This is the size of region in square degrees
            area = hp.nside2pixarea(nside=NSIDE, degrees=True) 
            if zmax is None:
                zmax = 1.0

        naming_dict = dict(ra='raJ2000',
                           dec='decJ2000',
                           redshift_true='redshift',
                           position_angle_true='morphology/positionAngle',
                           size_disk_true='morphology/diskHalfLightRadiusArcsec',
                           size_bulge_true='morphology/spheroidHalfLightRadiusArcsec',
                           size_minor_disk_true='morphology/diskMinorAxisArcsec',
                           size_minor_bulge_true='morphology/spheroidMinorAxisArcsec',
                           stellar_mass='totalMassStellar',
                           stellar_mass_disk='diskMassStellar')
        reverse_naming_dict = dict(list((naming_dict[key], key) for key in naming_dict))
        # print(reverse_naming_dict)
        # Veto galaxies used as main survey hosts as mddf hosts
        if veto_ms_gals:
            veto_galsdf = pd.read_csv(healpixelSN_fname)
            print(veto_galsdf.columns)
            veto_galsdf.rename(columns=dict(raJ2000_gal='ra', decJ2000_gal='dec'), inplace=True)
            print(veto_galsdf.columns)
            querystring = "ra > 52.479 and ra < 53.771 and dec < -27.533 and dec > -28.667"
            vetoed_galids = veto_galsdf.query(querystring).galaxy_id.values

        galsdf = pd.read_hdf(fname)
        if survey.lower() == 'mddf':
            galsdf = galsdf.query(querystring)

        galsdf.rename(columns=naming_dict, inplace=True)
  
        zdist = snsims.PowerLawRates(rng=np.random.RandomState(1 + randomSeedOffset),
                                     alpha=2.6e-5, beta=1.5, 
                                     fieldArea=area,
                                     surveyDuration=10.,
                                     cosmo=dc2Cosmo,
                                     zbinEdges=np.arange(0.001, zmax, 0.02))

        snPop_rng = np.random.RandomState(2 + randomSeedOffset)
        snPop = snsims.GMM_SALT2Params(numSN=None, zSamples=zdist.zSamples,
                                       rng=snPop_rng,
                                       mjdmin=59580, cosmo=dc2Cosmo,
                                       surveyDuration=10.)
    
    
        assert snPop.paramSamples.snid.max() - snPop.paramSamples.snid.size == -1

        max_redshift = galsdf.redshift.max()
        sn = DC2SN(galsdf, snPop, zmax=zmax, rng=np.random.RandomState(0+ randomSeedOffset))
        main_survey_mapper, hosted_sn_params = sn.assignHosts(binwidth=0.02,)
        if len(hosted_sn_params) > 0:
            hostedSNParamsPos = sn.get_positions(hosted_sn_params,
                                    np.random.RandomState(3 + randomSeedOffset))

            hostedSNParamsPos.snid = list(survey +
                                         '_{0}_{1}'.format(healpixelId, ind)
                                         for ind in
                                         hostedSNParamsPos.snid.values)

            out_name = os.path.join(args.out_dir, 'sn_{0}_{1}.csv'.format(healpixelId, survey))
            assert not os.path.isfile(out_name)
            hostedSNParamsPos.to_csv(out_name)
