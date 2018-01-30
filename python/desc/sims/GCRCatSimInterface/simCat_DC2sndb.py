"""
A module for classes that read in the protoDC2 catalog and assign SN to hosts

"""
from __future__ import division, print_function, absolute_import
__all__ = [ 'HostingSN']
import numpy as np
import pandas as pd
import snsims
from snsims import SersicSamples
from lsst.sims.utils import samplePatchOnSphere


class DC2SN(object):
    """ Class to associate SN in the SN params table with hosts
    
    Parameters
    ----------
    galsdf : pandas dataFrame
        provides galaxy table
    snpops : instance of SNPopulation
    redshiftvar : str, defaults to 'redshift'
    surveywidth : float
        width of the survey (5 degree for DC2)

    """
    def __init__(self,
                 galsdf,
                 snpops,
                 redshiftvar='redshift',
                 rng=np.random.RandomState(23),
                 assignBinWidth=0.02,
                 zmax=0.96,
                 surveyCenter=(55.064, -29.783),
                 surveyWidth=5.0,
                 hostprob=0.9):
        self.arcsec2Deg = 1.0 / 3600.
        self.snparams = snpops.paramSamples
        self.rng = rng
        self.zmax = zmax
        self.hostprob = hostprob
        self.redshiftvar = 'redshift'
        querystr = '{0} < {1}'.format(redshiftvar, zmax)
        self.galsdf = galsdf.query(querystr)
        
        # Sample of SN split into hosted/hostless
        hstd, hstless = self.splitByHosting(self.snparams, rng=self.rng,
                                            hostprob=self.hostprob, zmax=self.zmax)
        self.hostedSN = hstd
        self._hostedParams = None
        self._hostlessSN = hstless
        self.recalculate = False
        self.delta = surveyWidth / 2.0 
        self.centralRA = surveyCenter[0]
        self.centralDec = surveyCenter[1]
        self.ss = SersicSamples(self.rng)
    
    @staticmethod
    def splitByHosting(snparams, rng, hostprob=0.9, zmax=0.96):
        """Split the SN sample by hosting type
        """
        #snparams= pops.paramSamples
        snparams['rand_host'] = rng.uniform(size=len(snparams))
        hosted = snparams.query('z < @zmax and rand_host < @hostprob')
        hostlessSN = snparams.query('z >= @zmax or rand_host >= @hostprob')
        return hosted, hostlessSN

    @property
    def hostless_snparams(self):
        """
        calculate the parameter table for hostless SN
        """
        calculate = False
        if self._hostlessSN is None:
            calculate = True
        elif 'snra' not in self._hostlessSN.columns:
            calculate = True

        if calculate:
            ra , dec = samplePatchOnSphere(self.centralRA, self.centralDec, self.delta,
                                           size=len(self._hostlessSN))
            self._hostlessSN['snra'] = ra
            self._hostlessSN['sndec'] = dec
        return self._hostlessSN
    
    def assignHosts(self, binwidth=0.02):
        """
        Find hosts in redshift bins of 0.02 and populate them weighting
        them by stellar mass of host galaxies.
        
        Parameters
        ----------
        binwidth : float, defaults to 0.02
            width of redshift bin
        """
        galsdf = self.galsdf
        self.hostedSN['zbin'] = self.hostedSN.z // binwidth
        self.hostedSN['zbin'] = self.hostedSN.zbin.astype(np.int)
        galsdf['zbin'] = galsdf.redshift // binwidth
        galsdf['zbin'] = galsdf.zbin.astype(np.int)
        
        
        
        syslist = []

        for zbin in self.hostedSN.zbin.unique():
            querystring = 'zbin==@zbin'
            hostedtmp = self.hostedSN.query(querystring)
            numSN = hostedtmp.z.size
            gtmp  = galsdf.query(querystring)
            tot = gtmp.totalMassStellar.sum()
            gids = np.random.choice(gtmp.reset_index().galaxy_id, size=numSN, replace=False,
                                    p=gtmp.totalMassStellar/tot)
            #print(len(hostedtmp), len(gids))
            syslist.append(pd.DataFrame(dict(snid=hostedtmp.reset_index().snid, galaxy_id=gids)))
        
        joiner = pd.concat(syslist).set_index('snid')
        cols = self.hostedSN.columns
        keepcols = list(col for col in cols if col != 'z')
        params = self.hostedSN[keepcols].set_index('snid').join(joiner).set_index('galaxy_id').join(galsdf.set_index('galaxy_id'), rsuffix='_gals')
        #params = self.hostedSN.reset_index().set_index('snid').join(joiner).reset_index().set_index('galaxy_id').join(self.galsdf, rsuffix='_gals')
        params.rename(columns=dict(raJ2000='raJ2000_gal', decJ2000='decJ2000_gal', redshift='z'), inplace=True)
        
        k = joiner.reset_index().set_index('galaxy_id')
        print(k.columns)
        xx = params.join(k)
        cols = params.columns
        dropcols = ['index', 'z', 'zbin', 'zbin_gals', 'rand_host'] #+ list(a for a in cols if 'zbin' not in a)
        keepcols = list(col for col in cols if col not in dropcols)
        #xx = joiner.reset_index().set_index('galaxy_id').join(params)
        return joiner, xx
    
    def columndict(self):
        return dict(x0='Tx0', x1='Tx1', c='Tc', snra='raJ2000', sndec='decJ2000')

    def write(self, fname, paramsdf, redshiftvar='z'):
        x = paramsdf[redshiftvar].values 
        paramsdf.rename(columns=self.columndict, inplace=True)
        paramsdf['redshift'] = x
        paramsdf['Tredshift'] = x
        return
    
    def get_positions(self, paramsdf, rng, disk_a='size_disk_true', disk_b='morphology/diskMinorAxisArcsec',
                      bulge_a='size_disk_true', bulge_b='morphology/spheroidMinorAxisArcsec',
                      diskhalfLight='morphology/diskHalfLightRadiusArcsec',
                      bulgehalfLight='morphology/spheroidHalfLightRadiusArcsec',
                      positionAngle='morphology/positionAngle'):
        """
        Obtain ra dec of the SN and add it to the parameter table
        """
        pardf = paramsdf.copy()
        pardf['rand_SM'] = rng.uniform(size=len(pardf))
        pardf['diskmassratio'] = pardf['diskMassStellar'] / pardf['totalMassStellar']
        pardf['inDisk'] = pardf['diskmassratio'] > pardf['rand_SM']
        pardf.inDisk = pardf.inDisk.astype(np.int)
        x = pardf['inDisk']
        # Add angle samples
        pardf['diskAngle'] = self.ss.sampleAngles(pardf[disk_a], pardf[disk_b], numSamples=len(pardf))
        pardf['bulgeAngle'] = self.ss.sampleAngles(pardf[bulge_a], pardf[bulge_b], numSamples=len(pardf))
        
        pardf['diskRadius'] = self.ss.sampleRadius(pardf[diskhalfLight], numSamples=len(pardf), sersicIndex=1)
        pardf['bulgeRadius'] = self.ss.sampleRadius(pardf[bulgehalfLight], numSamples=len(pardf), sersicIndex=4)
        pardf['theta'] = pardf['diskAngle'] * x + (1-x) *pardf['bulgeAngle'] - pardf[positionAngle] +90 
        pardf['theta'] = np.radians(pardf.theta)
        pardf['radial'] = pardf['diskRadius'] * x + pardf['bulgeRadius']* (1-x)
        pardf['deltaRA'] = np.cos(pardf['theta']) * pardf['radial']
        pardf['deltaDec'] = np.sin(pardf['theta']) * pardf['radial']
        pardf['snra'] = pardf['raJ2000_gal'] + pardf['deltaRA'] * self.arcsec2Deg
        pardf['sndec'] = pardf['decJ2000_gal'] +pardf['deltaDec'] * self.arcsec2Deg
        return paramsdf.join(pardf[['snra', 'sndec']])
