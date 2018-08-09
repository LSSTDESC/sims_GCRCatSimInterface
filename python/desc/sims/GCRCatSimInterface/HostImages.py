import os
import numpy as np
import pandas as pd
from lsst.sims.utils import angularSeparation

__all__ = ['hostImage']

class hostImage(object):

    def __init__(self, ra_center, dec_center, fov):

        self.ra = ra_center
        self.dec = dec_center
        self.radius = fov


    def format_catalog(self, df_line, fits_file_name, image_dir):
        """Formats the output instance catalog to include entries for the FITS stamps produced by generate_lensed_host_***.py
		
		Parameters: 
		-----------
		df_line: string
			pandas data frame line
		fits_file_name: string
			the filename of the FITS stamp
		image_dir: string 
			the location of the FITS stamps 
        	"""
        fits_info = fits_file_name.split('_')
        lens_id = np.int(fits_info[0])
        sys_magNorm = np.float(fits_info[1])
        gal_type = fits_info[2].split('.')[0]

        galaxy_id = np.right_shift(lens_id, 10)
        galaxy_id *= 10000
        galaxy_id += 4*df_line['twinkles_system']
        sys_id = np.left_shift(galaxy_id, 10)
        if gal_type == 'bulge':
            sys_id += 97
        elif gal_type == 'disk':
            sys_id += 107
        else:
            raise ValueError('Invalid Galaxy Component Type in filename')

        cat_str = 'object %i %f %f %f %s %f 0 0 0 0 0 %s 0.01 0 CCM %f %f CCM %f %f\n' % (sys_id,
                                                                      df_line['raPhoSim_lens'],
                                                                      df_line['decPhoSim_lens'],
                                                                      sys_magNorm,
                                                                      df_line['sedFilepath'],
                                                                      df_line['redshift'],
                                                                      str(lens_id) + "_" + str(sys_magNorm) + "_" + str(gal_type) + ".fits",
                                                                      df_line['internalAv'],
                                                                      df_line['internalRv'],
                                                                      df_line['galacticAv'],
                                                                      df_line['galacticRv'])

       
        return cat_str

    def write_host_cat(self, image_dir, input_cat, output_cat, append=False):
        """Adds entries for each lensed host FITS stamp to output instance catalog
        Parameters:
        -----------
        image_dir: string
            the location of the FITS stamps 
        input_cat: string
            the list of FITS stamps to be read in
        output_cat: string
            the location of the output instance catalogs """
        
        host_df = pd.read_csv(input_cat)
        ang_sep_list = []
        image_list = os.listdir(image_dir)
        image_ids = np.array([image_name.split('_')[0] for image_name in image_list], dtype=np.int)
        
        for sys_ra, sys_dec in zip(host_df['raPhoSim_lens'], host_df['decPhoSim_lens']):
            ang_sep_list.append(angularSeparation(sys_ra, sys_dec, self.ra, self.dec))

        ang_sep = np.array(ang_sep_list)
        keep_idx = np.where(ang_sep < self.radius)
        host_image_df = host_df.iloc[keep_idx]

        unique_id_list = []

        if append:
            write_status = 'a'
        else:
            write_status = 'w'

        with open(output_cat, write_status) as f:
            for df_line in host_image_df.iterrows():
                line_uID = df_line[1]['uniqueId_lens']
                if line_uID in unique_id_list:
                    continue
                line_idx = np.where(image_ids == line_uID)[0]
                if len(line_idx) > 1:
                    raise ValueError('Multiple images have same unique lens id')
                line_name = image_list[line_idx[0]]
                f.write(self.format_catalog(df_line[1], line_name, image_dir))
                unique_id_list.append(line_uID)
