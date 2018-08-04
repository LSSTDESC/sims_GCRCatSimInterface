#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser(description='FITS stamp generator')
parser.add_argument('--outdir', type=str,
                    help='Output location for FITS stamps',
                    default=os.path.join(os.environ['SIMS_GCRCATSIMINTERFACE_DIR'], 
                                         'data', 'outputs'))

args = parser.parse_args()

from desc.sims.GCRCatSimInterface import generate_lensed_hosts_agn

generate_lensed_hosts_agn = generate_lensed_hosts_agn(args.outdir)

