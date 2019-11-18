import os
import sys
import glob
import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from get_config import get_config

try:
    config_file = sys.argv[1]
except IndexError:
    config_file = 'instcat_validation_config.ini'
config = get_config(config_file)

nullfmt = NullFormatter()

conn = sqlite3.connect(config['truthcat'])
df = pd.read_sql('select * from truth_summary', conn)
truth_cat_ids = set(df['id'])

run_description = config['run_description']
if run_description is None:
    run_description = ''

disable_MW = config['disable_MW']

#visits = dict(u=2338, g=159521, r=40325, i=479028, z=8005, y=5883)
visits = dict(i=6824)
bands = config['bands']

if ((config['frac_diff_min'] is not None) and
    (config['frac_diff_max'] is not None)):
    frac_diff_range = (config['frac_diff_min'], config['frac_diff_max'])
else:
    frac_diff_range = None


for i, band in enumerate(bands):
    visit = visits[band]
    print(visit, band)
    fig = plt.figure(figsize=(8, 6))
    visit_dir = f'{visit:08d}'
    flux_col = f'flux_{band}'
    tc_flux_col = flux_col
    if disable_MW:
        tc_flux_col += '_noMW'

    tc_flux_dict = dict([_ for _ in zip(df['id'], df[tc_flux_col])])

    tc_flux = []
    gal_flux = []
    ra = []
    dec = []
    gal_ids = set()

    my_truthcat_dir = config['ic_truthcat_prefix'] + f'_{visit_dir}'
    gal_cats = sorted(glob.glob(
        os.path.join(my_truthcat_dir, 'R*_S*_gal_truthcat.pickle')))

    for gal_cat in gal_cats:
        sys.stdout.write('.')
        sys.stdout.flush()
        gal_df = pd.read_pickle(gal_cat)
        gal_flux_dict = dict([_ for _ in zip(gal_df['gal_id'],
                                             gal_df[flux_col])])
        ra_dict = dict([_ for _ in zip(gal_df['gal_id'], gal_df['ra'])])
        dec_dict = dict([_ for _ in zip(gal_df['gal_id'], gal_df['dec'])])
        ids = truth_cat_ids.intersection(set(gal_df['gal_id']))
        for gal_id in ids:
            if gal_id not in gal_ids:
                tc_flux.append(tc_flux_dict[gal_id])
                gal_flux.append(gal_flux_dict[gal_id])
                ra.append(ra_dict[gal_id])
                dec.append(dec_dict[gal_id])
    print('!')
    tc_flux = np.array(tc_flux)
    gal_flux = np.array(gal_flux)*config['gal_flux_scale']
    frac_diff = (gal_flux - tc_flux)/tc_flux

    ra = np.array(ra)
    dec = np.array(dec)

    left, width = 0.1, 0.6
    bottom, height = 0.125, 0.8
    left_h = left + width + 0.02
    width_h = 0.2

    rect_hexbin = [left, bottom, width, height]
    rect_hist = [left_h, bottom, width_h, height]

    ax_hexbin = plt.axes(rect_hexbin)
    plt.xlabel(f'truthcat {flux_col} (Jy)')
    plt.ylabel('(instcat - truthcat)/truthcat')
    plt.title(f'visit {visit}, {flux_col} ' + run_description)
    log10_xmin = np.floor(np.log10(np.min(tc_flux)))
    log10_xmax = np.ceil(np.log10(np.max(tc_flux)))
    mappable = ax_hexbin.hexbin(tc_flux, frac_diff, xscale='log', gridsize=200,
                                mincnt=1, extent=(log10_xmin, log10_xmax,
                                                  *frac_diff_range))
    plt.colorbar(mappable)

    ax_hist = plt.axes(rect_hist)
    ax_hist.xaxis.set_major_formatter(nullfmt)
    ax_hist.hist(frac_diff, bins=200, orientation='horizontal',
                 range=frac_diff_range)

    plt.savefig(f'truthcat_vs_instcat_v{visit}-{band}_{run_description}.png')
