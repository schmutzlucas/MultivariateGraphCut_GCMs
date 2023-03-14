import os
import shutil
import xarray as xr
import numpy as np
from datetime import datetime

# Define the root directory containing all the model folders
root_dir = 'Y:\LucasSchmutz\MultivariateGraphCut_GCMs\download_day_unzip'

# Define the directory where the merged files will be stored
merged_dir = 'Y:\LucasSchmutz\MultivariateGraphCut_GCMs\download_day_merged'

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
                input_file_path = os.path.join(root_dir, model_dir, var_exp_dir, filename)
                input_files.append(input_file_path)

        # Print the input files being merged
        print(f"Merging {len(input_files)} files:")
        print('\n'.join(input_files))

        # Merge the netCDF files using xarray
        ds = xr.open_mfdataset(input_files, combine='nested', concat_dim='time')

        # Extract the start and end dates from the input files and convert them to datetime objects
        start_dates = []
        end_dates = []
        for input_file in input_files:
            time_obj = xr.open_dataset(input_file).time.values[0]
            start_date_str = np.datetime_as_string(time_obj, unit='D')
            start_date = datetime.strptime(start_date_str, '%Y-%m-%d')
            start_dates.append(start_date)

            time_obj = xr.open_dataset(input_file).time.values[-1]
            end_date_str = np.datetime_as_string(time_obj, unit='D')
            end_date = datetime.strptime(end_date_str, '%Y-%m-%d')
            end_dates.append(end_date)

        start_date = min(start_dates).strftime("%Y%m%d")
        end_date = max(end_dates).strftime("%Y%m%d")

        # Define output file name for merged file
        input_filename = os.path.basename(input_files[0])
        var = input_filename.split("_")[0]
        model = input_filename.split("_")[2]
        exp = input_filename.split("_")[3]
        ens = input_filename.split("_")[4]
        version = input_filename.split("_")[-1].split(".")[0]
        output_filename = f"{var}_{model}_{exp}_{ens}_{start_date}-{end_date}_{version}.nc"

        # Write merged file to network folder
        output_path = os.path.join(output_dir, output_filename)
        ds.to_netcdf(output_path)

# Remove the temporary directory if it exists
if os.path.exists(temp_dir):
    os.rmdir(temp_dir)
