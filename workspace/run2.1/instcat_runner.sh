out_dir="${@:1:1}"
config_file="${@:2:1}"
n_jobs="${@:3:1}"
echo 'out_dir '${out_dir}
echo 'config_file '${config_file}
for nn in "${@:3}";
do
    echo 'trying '${nn}
    python $SIMS_GCRCATSIMINTERFACE_DIR/bin.src/generateInstCat.py \
    --out_dir ${out_dir} \
    --config_file ${config_file} \
    --ids ${nn} --n_jobs ${n_jobs} --suppress_warnings \
    --host_data_dir $TWINKLES_DIR/data/ --enable_sprinkler &
done
wait
