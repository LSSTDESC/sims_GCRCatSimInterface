out_dir=$SCRATCH/instcat_test/
obs=2333

python ../../bin.src/generateInstCat.py \
--out_dir ${out_dir} \
--config_file config_file_2.1.wfd.json \
--ids ${obs} --n_jobs 1 --suppress_warnings \

wait
