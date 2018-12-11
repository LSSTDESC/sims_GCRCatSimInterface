sne_dir=$SCRATCH/cosmoDC2_v1.1.4_sne/
python write_healpixel.py --out_dir ${sne_dir}sne_csv_mddf \
--data_root ${sne_dir}gal_hdf5_mddf \
--survey mDDF --zmax 1.4 --healpixelId 9043
