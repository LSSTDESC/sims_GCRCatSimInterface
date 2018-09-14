#!/usr/bin/env python
import argparse
import warnings
import os
import copy
import time
import multiprocessing
import numbers
from astropy._erfa import ErfaWarning

with warnings.catch_warnings():
    warnings.filterwarnings('ignore', '\nThis call', UserWarning)
    warnings.filterwarnings('ignore', 'Duplicate object type', UserWarning)
    warnings.filterwarnings('ignore', 'No md5 sum', UserWarning)
    warnings.filterwarnings('ignore', 'ERFA function', ErfaWarning)
    warnings.filterwarnings('ignore', 'Numpy has detected', FutureWarning)
    warnings.filterwarnings('ignore', 'divide by zero', RuntimeWarning)
    warnings.filterwarnings('ignore', 'invalid value', RuntimeWarning)

    from desc.sims.GCRCatSimInterface import InstanceCatalogWriter


def generate_instance_catalog(args=None, lock=None):

    with warnings.catch_warnings():
        if args.suppress_warnings:
            warnings.filterwarnings('ignore', '\nThis call', UserWarning)
            warnings.filterwarnings('ignore', 'Duplicate object type', UserWarning)
            warnings.filterwarnings('ignore', 'No md5 sum', UserWarning)
            warnings.filterwarnings('ignore', 'ERFA function', ErfaWarning)
            warnings.filterwarnings('ignore', 'Numpy has detected', FutureWarning)
            warnings.filterwarnings('ignore', 'divide by zero', RuntimeWarning)
            warnings.filterwarnings('ignore', 'invalid value', RuntimeWarning)

        if not hasattr(generate_instance_catalog, 'instcat_writer'):
            instcat_writer = InstanceCatalogWriter(args.db, args.descqa_catalog,
                                                   dither=not args.disable_dithering,
                                                   min_mag=args.min_mag,
                                                   minsource=args.minsource,
                                                   proper_motion=args.enable_proper_motion,
                                                   protoDC2_ra=args.protoDC2_ra,
                                                   protoDC2_dec=args.protoDC2_dec,
                                                   agn_db_name=args.agn_db_name,
                                                   sn_db_name=args.sn_db_name,
					           host_image_dir=args.host_image_dir,
                                                   host_data_dir=args.host_data_dir,
                                                   sprinkler=args.enable_sprinkler)

            generate_instance_catalog.instcat_writer = instcat_writer


        for obsHistID in args.ids:
            if args.job_log is not None:
                if lock is not None:
                    lock.acquire()
                with open(args.job_log, 'a') as out_file:
                    out_file.write('starting %d at time %.0f\n' % (obsHistID, time.time()))
                if lock is not None:
                    lock.release()

            full_out_dir = os.path.join(args.out_dir, '%.8d' % obsHistID)

            generate_instance_catalog.instcat_writer.write_catalog(obsHistID,
                                                                   out_dir=full_out_dir,
                                                                   fov=args.fov,
                                                                   status_dir=args.out_dir)

            if args.job_log is not None:
                if lock is not None:
                    lock.acquire()
                with open(args.job_log, 'a') as out_file:
                    out_file.write('ending %d at time %.0f\n' % (obsHistID, time.time()))
                if lock is not None:
                    lock.release()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Instance catalog generator')
    parser.add_argument('--db', type=str,
                        help='path to the OpSim database to query')
    parser.add_argument('--agn_db_name', type=str,
                        help='File of AGN parameters generated by create_agn_db.py')
    parser.add_argument('--sn_db_name', type=str, default=None,
                        help='File of SN parameters generated by create_sne_db.py')
    parser.add_argument('--host_image_dir', type=str,
                        help='Location of FITS stamps of lensed host images produced by generate_lensed_hosts_***.py',
                        default=os.path.join(os.environ['SIMS_GCRCATSIMINTERFACE_DIR'], 'data', 'outputs'))
    parser.add_argument('--host_data_dir', type=str,
                        help='Name of csv file of lensed host data created by the sprinkler.',
                        default=os.path.join(os.environ['SIMS_GCRCATSIMINTERFACE_DIR'],'data'))
    parser.add_argument('--descqa_catalog', type=str, default='protoDC2',
                        help='the desired DESCQA catalog')
    parser.add_argument('--out_dir', type=str, help='directory where output will be written',
                        default=os.path.join(os.environ['SIMS_GCRCATSIMINTERFACE_DIR'],'data', 'outputs'))
    parser.add_argument('--ids', type=int, nargs='+',
                        default=None,
                        help='obsHistID to generate InstanceCatalog for (a list)')
    parser.add_argument('--disable_dithering', default=False,
                        action='store_true',
                        help='flag to disable dithering')
    parser.add_argument('--min_mag', type=float, default=10.0,
                        help='the minimum magintude for stars')
    parser.add_argument('--fov', type=float, default=2.0,
                        help='field of view radius in degrees')
    parser.add_argument('--enable_proper_motion', default=False,
                        action='store_true',
                        help='flag to enable proper motion')
    parser.add_argument('--minsource', type=int, default=100,
                        help='mininum #objects in a trimmed instance catalog')
    parser.add_argument('--protoDC2_ra', type=float, default=0,
                        help='RA (J2000 degrees) of the new protoDC2 center')
    parser.add_argument('--protoDC2_dec', type=float, default=0,
                        help='Dec (J2000 degrees) of the new protoDC2 center')
    parser.add_argument('--enable_sprinkler', default=False, action='store_true',
                        help='flag to enable the sprinkler')
    parser.add_argument('--suppress_warnings', default=False, action='store_true',
                        help='flag to suppress warnings')
    parser.add_argument('--n_jobs', type=int, default=1,
                        help='Number of jobs to run in parallel with multiprocessing')
    parser.add_argument('--job_log', type=str, default=None,
                        help="file where we will write 'job started/completed' messages")
    args = parser.parse_args()

    if args.n_jobs==1 or isinstance(args.ids, numbers.Number) or len(args.ids)==1:
        generate_instance_catalog(args=args)
    else:
        lock = multiprocessing.Lock()
        job_list = []
        n_id = len(args.ids)//args.n_jobs  # number of ids per job
        for i_start in range(0, len(args.ids), n_id):
            local_args = copy.deepcopy(args)
            local_args.ids = args.ids[i_start:i_start+n_id]
            p = multiprocessing.Process(target=generate_instance_catalog,
                                        kwargs={'args':local_args, 'lock':lock})
            p.start()
            job_list.append(p)

        for p in job_list:
            p.join()

    with open(args.job_log, 'a') as out_file:
        out_file.write('%s should be completed\n' % str(args.ids))
