import os
import sys
import glob
import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

flux_col = 'flux_i'

conn = sqlite3.connect('truth_cat/summary_table_hp9430_n100000.sqlite3')
df = pd.read_sql('select * from truth_summary', conn)
truth_cat_ids = set(df['id'])
tc_flux_dict = dict([_ for _ in zip(df['id'], df[flux_col])])

tc_flux = []
gal_flux = []
ra = []
dec = []
gal_ids = set()

gal_cats = sorted(glob.glob('my_truth_cats/R*_S*_gal_truthcat.pickle'))
for gal_cat in gal_cats:
    sys.stdout.write('.')
    gal_df = pd.read_pickle(gal_cat)
    gal_flux_dict = dict([_ for _ in zip(gal_df['gal_id'], gal_df[flux_col])])
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
gal_flux = np.array(gal_flux)
frac_diff = (gal_flux - tc_flux)/tc_flux

ra = np.array(ra)
dec = np.array(dec)

nullfmt = NullFormatter()

left, width = 0.1, 0.6
bottom, height = 0.125, 0.8
left_h = left + width + 0.02
width_h = 0.2

rect_hexbin = [left, bottom, width, height]
rect_hist = [left_h, bottom, width_h, height]

plt.figure(figsize=(7, 4))

ax_hexbin = plt.axes(rect_hexbin)
plt.xlabel('truthcat flux_i (Jy)')
plt.ylabel('(instcat - truthcat)/truthcat')
plt.title('visit 479028, flux_i')
mappable = ax_hexbin.hexbin(tc_flux, frac_diff, xscale='log', gridsize=200,
                            mincnt=1)
plt.colorbar(mappable)

ax_hist = plt.axes(rect_hist)
ax_hist.xaxis.set_major_formatter(nullfmt)
ax_hist.hist(frac_diff, bins=200, orientation='horizontal')
