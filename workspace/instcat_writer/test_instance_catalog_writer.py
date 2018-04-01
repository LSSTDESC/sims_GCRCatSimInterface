from desc.sims.GCRCatSimInterface import InstanceCatalogWriter
import os

if __name__ == "__main__":

    opsim_db_name = os.path.join(os.environ['SCRATCH'],'OpSimData',
                                 'minion_1016_sqlite_new_dithers.db')

    assert os.path.exists(opsim_db_name)

    agn_db_name = os.path.join(os.environ['SCRATCH'], 'proto_dc2_agn',
                               'test_agn.db')

    assert os.path.exists(agn_db_name)

    ic_writer = InstanceCatalogWriter(opsim_db_name, 'proto-dc2_v3.0',
                                      protoDC2_ra=53.0, protoDC2_dec=-28.0,
                                      agn_db_name=agn_db_name,
                                      sprinkler=True,
                                      dither=False)

    out_dir = 'sprinkled'
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    assert os.path.isdir(out_dir)

    ic_writer.write_catalog(230, out_dir=out_dir, fov=0.1)
