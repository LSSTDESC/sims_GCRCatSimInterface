from ExtractSNeParams import add_unsprinkled_sne_params
import os

db_dir = '/astro/store/pogo3/danielsf/truth_181002'
assert os.path.isdir(db_dir)

db_file = os.path.join(db_dir, 'truth_params_181002_protodc2_v3.db')
assert os.path.isfile(db_file)

add_unsprinkled_sne_params(db_file)
