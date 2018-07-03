import os
from desc.sims.GCRCatSimInterface import write_sprinkled_truth_db
from lsst.sims.utils import ObservationMetaData

agn_dir = '/global/projecta/projectdirs/lsst/groups/SSim/DC2'
assert os.path.isdir(agn_dir)

agn_db = os.path.join(agn_dir, 'agn_db_mbh_7.0_m_i_30.0.sqlite')
assert os.path.isfile(agn_db)

obs = ObservationMetaData(pointingRA=53.13231600394978216,
                          pointingDec=-28.03064232384672749,
                          mjd=62746.27986361111106817,
                          boundType='circle',
                          boundLength=0.5)

write_sprinkled_truth_db(obs, agn_db=agn_db)
