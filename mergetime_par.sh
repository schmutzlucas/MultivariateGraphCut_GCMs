#!/bin/bash

# Define the root directory containing all the model folders
root_dir="Y:/scratch_idyst/LucasSchmutz/MultivariateGraphCut_GCMs/download_day_unzip"

# Define the directory where the merged files will be stored
merged_dir="Y:/scratch_idyst/LucasSchmutz/MultivariateGraphCut_GCMs/download_day_merged"

# Define the number of parallel processes to use
nprocs=4

# Loop over all model folders
for model_dir in "$root_dir"/*; do
    # Ignore any non-directory files in the root directory
    if [ ! -d "$model_dir" ]; then
        continue
    fi

    # Loop over all variable/experiment subfolders in the model directory
    for var_exp_dir in "$model_dir"/*; do
        # Ignore any non-directory files in the variable/experiment subfolder
        if [ ! -d "$var_exp_dir" ]; then
            continue
        fi

        # Construct the output directory path for the merged files
        output_dir="$merged_dir/$(basename "$model_dir")/$(basename "$var_exp_dir")"
        mkdir -p "$output_dir"

        # Loop over all netCDF files in the variable/experiment subfolder
        input_files=()
        for filename in "$var_exp_dir"/*.nc; do
            if [ -f "$filename" ]; then
                input_files+=("$filename")
            fi
        done

        # Merge the netCDF files using cdo mergetime in parallel
        output_filename="$(basename "${input_files[0]}")"
        output_filename="${output_filename/_19500101-/}"
        output_filename="${output_filename/_v/_merged_v}"
        output_path="$output_dir/$output_filename"
        echo "${input_files[@]}" | tr ' ' '\n' | \
            parallel -j $nprocs --no-notice cdo mergetime {} "$output_path"
    done
done
