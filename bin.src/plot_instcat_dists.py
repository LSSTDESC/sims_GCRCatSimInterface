#!/usr/bin/env python
"""
Script to plot distributions of various quantities in a set of
phosim instance catalog files produced by generateInstCat.py.
"""
import os
import gzip
import argparse
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

class HistArray(object):
    "Class to manage subplotting of histograms."
    def __init__(self, title='', figsize=(9, 12), shape=(5, 2)):
        """
        Parameters
        ----------
        title: str ['']
            Overall title of the figure.
        figsize: tuple [(9, 12)]
            Figure dimensions in x, y inches
        shape: tuple [(5, 2)]
            Number of subplots in the y and x dimensions, respectively.
        """
        plt.rcParams['figure.figsize'] = figsize
        self.fig = plt.figure()
        frame_axes = self.fig.add_subplot(111, frameon=False, xticklabels=(),
                                          yticklabels=())
        frame_axes.set_title(title)
        self.shape = shape
        self.num_subplots = 0

    def plot_hists(self, columns, xlabel, yscale='log', bins=50):
        """
        Plot histograms for a given columnar quantity, e.g., RA, Dec, etc..

        Parameters
        ----------
        columns: dict of lists
            Dictionary of lists of column values, keyed by object file name
            (e.g., 'gal_cat_197356').
        xlabel: str
            Label of x-axis.
        yscale: str ['log']
            y-axis scaling, 'log' or 'lin'
        bins: int [50]
            Number of bins for histogram.
        """
        self.num_subplots += 1
        self.fig.add_subplot(self.shape[0], self.shape[1], self.num_subplots)
        column_values = np.concatenate(tuple(columns.values()))
        x_range = min(column_values), max(column_values)
        for item in columns:
            plt.hist(columns[item], range=x_range, bins=bins, histtype='step',
                     label=item.split('.')[0])
        plt.xlabel(xlabel)
        plt.yscale(yscale)
        plt.legend(loc=1, fontsize=6)

def plot_instcat_dists(phosim_file, figsize=(9, 12)):
    """
    Create a multipanel plot of histograms of various columns in the
    phosim instance catalog.

    Parameters
    ----------
    phosim_file: str
        Instance catalog file containing includeobj references for each
        object type.
    figsize: tuple [(9, 12)]
        Figure dimensions in x, y inches
    """
    instcat_dir = os.path.split(phosim_file)[0]
    hist_array = HistArray(title=phosim_file, figsize=figsize)
    object_files = []
    with open(phosim_file) as phosim_input:
        for line in phosim_input:
            if line.startswith('includeobj'):
                object_files.append(line.strip().split()[-1])
    ra = defaultdict(list)
    dec = defaultdict(list)
    magnorm = defaultdict(list)
    major_axis = defaultdict(list)
    minor_axis = defaultdict(list)
    num_zero_major = 0
    num_zero_minor = 0
    axis_ratio = defaultdict(list)
    redshift = defaultdict(list)
    gamma1 = defaultdict(list)
    gamma2 = defaultdict(list)
    kappa = defaultdict(list)
    for item in object_files:
        with gzip.open(os.path.join(instcat_dir, item), 'r') as objects:
            for line in objects:
                tokens = line.split()
                ra[item].append(float(tokens[2]))
                dec[item].append(float(tokens[3]))
                if float(tokens[4]) < 1000:
                    magnorm[item].append(float(tokens[4]))
                if float(tokens[6]) > 0:
                    redshift[item].append(float(tokens[6]))
                if 'sersic2d' in str(line):
                    major_axis[item].append(float(tokens[13]))
                    minor_axis[item].append(float(tokens[14]))
                    if major_axis[item][-1] <= 0:
                        num_zero_major += 1
                    if minor_axis[item][-1] <= 0:
                        num_zero_minor += 1
                    else:
                        axis_ratio[item].append(major_axis[item][-1]/
                                                minor_axis[item][-1])
                        if axis_ratio[item][-1] > 1000:
                            print(line.strip())
                    gamma1[item].append(float(tokens[7]))
                    gamma2[item].append(float(tokens[8]))
                    kappa[item].append(float(tokens[9]))
    hist_array.plot_hists(ra, xlabel='RA (degrees)')
    hist_array.plot_hists(dec, xlabel='Dec (degrees)')
    hist_array.plot_hists(magnorm, xlabel='magnorm')
    hist_array.plot_hists(redshift, xlabel='redshift')
    hist_array.plot_hists(major_axis, xlabel='sersic2d major axis (arcsec)')
    hist_array.plot_hists(minor_axis, xlabel='sersic2d minor axis (arcsec)')
    hist_array.plot_hists(axis_ratio, xlabel='major/minor (#<=0:{} major, {} minor)'.format(num_zero_major, num_zero_minor))
    hist_array.plot_hists(gamma1, xlabel='gamma1')
    hist_array.plot_hists(gamma2, xlabel='gamma2')
    hist_array.plot_hists(kappa, xlabel='kappa')
    plt.tight_layout()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="plot distributions of quantities from phosim instance catalogs")
    parser.add_argument('phosim_file', type=str,
                        help='Instance catalog file with includeobj references')
    parser.add_argument('--pngfile', type=str, default=None,
                        help='Output png file. If None, then the name will be '
                        'derived from the instance catalog file name')
    args = parser.parse_args()

    plot_instcat_dists(args.phosim_file)
    pngfile = args.pngfile
    if pngfile is None:
        pngfile = os.path.basename(args.phosim_file).replace('.txt', '.png')
    plt.savefig(pngfile)
