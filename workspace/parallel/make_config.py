import json
import os

project_dir = '/global/projecta/projectdirs/lsst/groups/SSim/DC2'

config_dict = {}

config_dict['db'] = os.path.join(project_dir,'minion_1016_desc_dithered_v4.db')
config_dict['agn_db_name'] =  os.path.join(project_dir, 'agn_db_mbh_7.0_m_i_30.0_cosmodc2_180912.db')
config_dict['sn_db_name'] =  os.path.join(project_dir, 'sne_params_wfd_cosmodc2_180911.db')
config_dict['descqa_catalog'] = 'cosmoDC2_v1.0_image_addon_knots'
config_dict['fov'] = 0.02
config_dict['n_jobs'] = 2

with open('config_file_test.json', 'w') as out_file:
    json.dump(config_dict, out_file)
