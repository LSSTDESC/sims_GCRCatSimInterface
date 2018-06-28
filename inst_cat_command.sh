output_dir=workspace/catalogs/truth_test/

python bin.src/generateInstCat.py --db /global/projecta/projectdirs/lsst/groups/SSim/DC2/minion_1016_desc_dithered_v4.db --agn_db_name /global/projecta/projectdirs/lsst/groups/SSim/DC2/agn_db_mbh_7.0_m_i_30.0.sqlite --descqa_catalog proto-dc2_v4.6.1 --fov 0.5 --enable_sprinkler --out_dir ${output_dir} --protoDC2_ra 55.064 --protoDC2_dec -29.783 --ids 2145656
