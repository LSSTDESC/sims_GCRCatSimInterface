import argparse
import os
from lsst.sims.utils import ObservationMetaData
from lsst.sims.photUtils import BandpassDict
from desc.sims.GCRCatSimInterface import write_sprinkled_lc


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--param_file', type=str, default=None,
                        help='Path to database created by '
                        'generate_truth_params.py')
    parser.add_argument('--ptng_dir', type=str, default=None,
                        help='Directory containing text files which list '
                        'the obsHistIDs of the pointings being simulated')
    parser.add_argument('--opsim_db', type=str, default=None,
                        help='Path to the OpSim database of pointings being '
                        'simulated')
    parser.add_argument('--out_file', type=str, default=None,
                        help='Name of light curve database to be written '
                        '(must not already exist)')
    parser.add_argument('--RA', type=float, default=55.064,
                        help='RA at center of simulated area in degrees '
                        '(default 55.064)')
    parser.add_argument('--Dec', type=float, default=-29.783,
                        help='Dec at center of simulated area in degrees '
                        '(default -29.783)')
    parser.add_argument('--fov', type=float, default=4.0,
                        help='radius in degrees of simulated area '
                        '(default 4.0)')

    args = parser.parse_args()
    assert os.path.isfile(args.param_file)
    assert os.path.isdir(args.ptng_dir)
    assert os.path.isfile(args.opsim_db)
    assert not os.path.isfile(args.out_file)

    obs_tot = ObservationMetaData(pointingRA=args.RA, pointingDec=args.Dec,
                                  boundType='circle', boundLength=args.fov)

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

    write_sprinkled_lc(args.out_file, obs_tot,
                       args.ptng_dir, args.opsim_db,
                       sql_file_name=args.param_file,
                       bp_dict=bp_dict)
