"""
"""
__all__ = ['DC2SN', 'DC2']
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
# from varpop import SersicSamples
from snsims import SersicSamples
from lsst.sims.utils import samplePatchOnSphere
# from opsimsummary import Tiling

DC2 = FlatLambdaCDM(H0=71, Om0=0.265, Ob0=0.0448)


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
                                            hostprob=self.hostprob,
                                            zmax=self.zmax)
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
    
    def assignHosts(self, binwidth=0.02,
                    ra_center=None,
                    dec_center=None,
                    ra_dec_width=None):
        """
        Find hosts in redshift bins of 0.02 and populate them weighting
        them by stellar mass of host galaxies.
        
        Parameters
        ----------
        binwidth : float, defaults to 0.02
            width of redshift bin

        ra_center, dec_center, ra_dec_width are in degrees;
        denote the center and size of the DDF
        """
        galsdf = self.galsdf
        self.hostedSN['zbin'] = self.hostedSN.z // binwidth
        self.hostedSN['zbin'] = self.hostedSN.zbin.astype(np.int)
        galsdf['zbin'] = galsdf.redshift // binwidth
        galsdf['zbin'] = galsdf.zbin.astype(np.int)

        syslist = []

        for zbin in self.hostedSN.zbin.unique():
            print("This is zbin", zbin)
            querystring = 'zbin==@zbin'
            hostedtmp = self.hostedSN.query(querystring)
            numSN = hostedtmp.z.size
            gtmp  = galsdf.query(querystring)
            if len(gtmp) == 0:
                continue
            p = gtmp.totalMassStellar.values
            p /= p.sum()
            print('the total is ', p.sum())

            print('numSN %d len(gtmp) %d' % (numSN, len(gtmp)))
            if numSN>len(gtmp):
                raise RuntimeError('numSN %d len(gtmp) %d' %
                                   (numSN, len(gtmp)))

            g_candidates = gtmp.reset_index()
            gid_candidates = g_candidates.galaxy_id.values
            ra_candidates = g_candidates.raJ2000.values
            dec_candidates = g_candidates.decJ2000.values

            gid_idx = self.rng.choice(np.arange(len(gid_candidates), dtype=int),
                                      size=numSN, replace=False,
                                      p=p) #gtmp.totalMassStellar/tot)

            gids = gid_candidates[gid_idx]
            ra_candidates = ra_candidates[gid_idx]
            dec_candidates = dec_candidates[gid_idx]

            if ra_center is not None:
                ra_min = ra_center-ra_dec_width
                ra_max = ra_center+ra_dec_width
                dec_min = dec_center-ra_dec_width
                dec_max = dec_center+ra_dec_width
                # map ra onto an actually rectangular coordinate
                # so that we can keep all of the points in the DDF
                ra_candidates = (ra_center
                  +(ra_candidates-ra_center)*np.cos(np.radians(dec_candidates)))

                valid = np.where(np.logical_and(ra_candidates >= ra_min,
                                 np.logical_and(ra_candidates <= ra_max,
                                 np.logical_and(dec_candidates >= dec_min,
                                                dec_candidates <= dec_max))))

                gids = gids[valid]
                hostedtmp = hostedtmp.iloc[valid]

            print(len(hostedtmp), len(gids))
            syslist.append(pd.DataFrame(dict(snid=hostedtmp.reset_index().snid, galaxy_id=gids)))

        if len(syslist) == 0:
            return [], []
        joiner = pd.concat(syslist).set_index('snid')
        cols = self.hostedSN.columns
        keepcols = list(col for col in cols if col != 'z')
        params = self.hostedSN[keepcols].set_index('snid').join(joiner).set_index('galaxy_id').join(galsdf.set_index('galaxy_id'), rsuffix='_gals')
        #params = self.hostedSN.reset_index().set_index('snid').join(joiner).reset_index().set_index('galaxy_id').join(self.galsdf, rsuffix='_gals')
        params.rename(columns=dict(raJ2000='raJ2000_gal', decJ2000='decJ2000_gal', redshift='z'), inplace=True)
        
        k = joiner.reset_index().set_index('galaxy_id')
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

        # Is the SN in the disk or bulge ?
        # We will have a column `inDisk` which is 1 if the SN is in the disk, 0 otherwise.
        pardf['rand_SM'] = rng.uniform(size=len(pardf))
        pardf['diskmassratio'] = pardf['diskMassStellar'] / pardf['totalMassStellar']
        pardf['inDisk'] = pardf['diskmassratio'] > pardf['rand_SM']
        pardf.inDisk = pardf.inDisk.astype(np.int)
        x = pardf['inDisk']
        # Add angle samples
        # pardf['diskAngle'] = self.ss.sampleAngles(pardf[disk_a], pardf[disk_b], numSamples=len(pardf))
        pardf['diskAngle'] = self.ss.sampleAngles(pardf[diskhalfLight], pardf[disk_b], numSamples=len(pardf))
        pardf['bulgeAngle'] = self.ss.sampleAngles(pardf[bulgehalfLight], pardf[bulge_b], numSamples=len(pardf))
        # pardf['bulgeAngle'] = self.ss.sampleAngles(pardf[bulge_a], pardf[bulge_b], numSamples=len(pardf))
        
        pardf['diskRadius'] = self.ss.sampleRadius(pardf[diskhalfLight], numSamples=len(pardf), sersicIndex=1)
        pardf['bulgeRadius'] = self.ss.sampleRadius(pardf[bulgehalfLight], numSamples=len(pardf), sersicIndex=4)
        pardf['theta'] = pardf['diskAngle'] * x + (1-x) *pardf['bulgeAngle'] - pardf[positionAngle] + 90 
        pardf['theta'] = np.radians(pardf.theta)
        pardf['radial'] = pardf['diskRadius'] * x + pardf['bulgeRadius']* (1-x)
        pardf['deltaRA'] = np.cos(pardf['theta']) * pardf['radial']
        pardf['deltaDec'] = np.sin(pardf['theta']) * pardf['radial']
        pardf['snra'] = pardf['raJ2000_gal'] + pardf['deltaRA'] * self.arcsec2Deg
        pardf['sndec'] = pardf['decJ2000_gal'] +pardf['deltaDec'] * self.arcsec2Deg
        return paramsdf.join(pardf[['snra', 'sndec']])
