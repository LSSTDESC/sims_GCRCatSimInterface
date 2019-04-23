out_dir="${@:1:1}"
config_file="${@:2:1}"
n_parallel_jobs="${@:3:1}"
obsid_list=("${@:4}")
n_obsid=${#obsid_list[@]}

echo 'out_dir '${out_dir}
echo 'config_file '${config_file}
declare -i n_per_job
n_per_job=${n_obsid}/${n_parallel_jobs}

echo 'n_per_job '${n_per_job}

declare -i i_start
i_start=0

while [ $i_start -lt $n_obsid ];
do
    id_list=${obsid_list[@]:${i_start}:${n_per_job}}
    i_start+=${n_per_job}

    python -W'ignore' $SIMS_GCRCATSIMINTERFACE_DIR/bin.src/generateInstCat.py \
    --out_dir ${out_dir} \
    --config_file ${config_file} \
    --ids ${id_list[@]} --n_jobs ${n_parallel_jobs} --suppress_warnings \
    --host_data_dir $TWINKLES_DIR/data/ --enable_sprinkler &
done
wait
