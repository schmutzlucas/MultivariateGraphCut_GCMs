import os
import shutil
import xarray as xr
import numpy as np
from datetime import datetime

# Define the root directory containing all the model folders
root_dir = 'Y:\LucasSchmutz\MultivariateGraphCut_GCMs\download_day_unzip'

# Define the directory where the merged files will be stored
merged_dir = 'Y:\LucasSchmutz\MultivariateGraphCut_GCMs\download_day_merged_multithread'

# Loop over all model folders
for model_dir in os.listdir(root_dir):
    # Ignore any non-directory files in the root directory
    if not os.path.isdir(os.path.join(root_dir, model_dir)):
        continue

    # Loop over all variable/experiment subfolders in the model directory
    for var_exp_dir in os.listdir(os.path.join(root_dir, model_dir)):
        # Ignore any non-directory files in the variable/experiment subfolder
        if not os.path.isdir(os.path.join(root_dir, model_dir, var_exp_dir)):
            continue

        # Construct the output directory path for the merged files
        output_dir = os.path.join(merged_dir, model_dir, var_exp_dir)
        os.makedirs(output_dir, exist_ok=True)

        # Loop over all netCDF files in the variable/experiment subfolder
        input_files = []
        for filename in os.listdir(
                os.path.join(root_dir, model_dir, var_exp_dir)):
            if filename.endswith(".nc"):
                input_file_path = os.path.join(root_dir, model_dir, var_exp_dir,
                                               filename)
                input_files.append(input_file_path)

        # Print the input files being merged
        print(f"Merging {len(input_files)} files:")
        print('\n'.join(input_files))
        # Get the current system time
        now = datetime.now()

        # Print the current system time
        print("Current time:", now)

        # Merge the netCDF files using xarray
        try:
            ds = xr.open_mfdataset(input_files, combine='nested',
                                   concat_dim='time')

            # Extract the start and end dates from the input files and convert
            # them to datetime objects
            start_date = np.datetime_as_string(ds.time.values[0],
                                               unit='D').replace('-', '')
            end_date = np.datetime_as_string(ds.time.values[-1],
                                             unit='D').replace('-', '')

            # Define output file name for merged file
            input_filename = os.path.basename(input_files[0])
            var = input_filename.split("_")[0]
            model = input_filename.split("_")[2]
            exp = input_filename.split("_")[3]
            ens = input_filename.split("_")[4]
            version = input_filename.split("_")[-1].split(".")[0]
            output_filename = f"{var}_{model}_{exp}_{ens}_{start_date}-{end_date}_{version}.nc"

            # Check if output file already exists before writing
            output_path = os.path.join(output_dir, output_filename)
            print(output_path)
            if os.path.exists(output_path):
                print(
                    f"Output file {output_filename} already exists, skipping...")
            else:
                # Write merged file to network folder
                ds.to_netcdf(output_path)
        except Exception as e:
            print(f"An error occurred while processing {input_files}: {e}")

# Remove the temporary directory if it exists
if os.path.exists(temp_dir):
    os.rmdir(temp_dir)
