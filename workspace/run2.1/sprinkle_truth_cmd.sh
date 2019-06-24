#for obs in 277065 159489 40320 184884 1473 7372;
#for obs in 538392;
#for obs in 1405954;
#for obs in 196570;
#for obs in 68022 219928 1037341 2141478;
#for obs in 306189;

out_dir=$SCRATCH/sprinkled_truth_example_190624/

for obs in 2500100;
do
    /usr/bin/time --verbose python ../../bin.src/generateInstCat.py \
    --out_dir ${out_dir} \
    --config_file config_sprinkled_truth.json \
    --ids ${obs} --n_jobs 1 --suppress_warnings \
    --enable_sprinkler --host_data_dir $TWINKLES_DIR/data/ 

done

wait
