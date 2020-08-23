sne_dir=$SCRATCH/cosmoDC2_v1.1.4_sne/
python write_healpixel.py --out_dir ${sne_dir}sne_csv --data_root ${sne_dir}gal_hdf5 --survey MS --zmax 1.0
