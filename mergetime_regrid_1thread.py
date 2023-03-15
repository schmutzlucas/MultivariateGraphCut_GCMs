import os
import shutil
import xarray as xr
import numpy as np
from datetime import datetime
import xesmf as xe

# Define the root directory containing all the model folders
root_dir = 'Y:\LucasSchmutz\MultivariateGraphCut_GCMs\download_day_unzip'

# Define the directory where the merged files will be stored
merged_dir = 'Y:\LucasSchmutz\MultivariateGraphCut_GCMs\download_day_merged_regrid_1thread'

# Define the path to the regridding grid file
grid_path = 'my_grid.txt'

try:
    # Load the grid information
    with open(grid_path) as f:
        lines = f.readlines()
        xsize = int(lines[1].split('=')[1].strip())
        ysize = int(lines[2].split('=')[1].strip())
        xfirst = int(lines[3].split('=')[1].strip())
        xinc = int(lines[4].split('=')[1].strip())
        yfirst = int(lines[5].split('=')[1].strip())
        yinc = int(lines[6].split('=')[1].strip())
except Exception as e:
    print(f"Error occurred while loading grid information: {e}")
    exit(1)

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
        print(f"Merging and regridding {len(input_files)} files:")
        print('\n'.join(input_files))
        # Get the current system time
        now = datetime.now()

        # Print the current system time
        print("Current time:", now)

        try:
            # Merge the netCDF files using xarray
            ds = xr.open_mfdataset(input_files, combine='nested', concat_dim='time')

            # Regrid the dataset to the destination grid
            regridder = xe.Regridder(ds, ds_out, 'bilinear')
            ds_regrid = regridder(ds)

            # Extract the start and end dates from the input files and convert
            # them to datetime objects
            start_date = np.datetime_as_string(ds_regrid.time.values[0],
                                               unit='D').replace('-', '')
            end_date = np.datetime_as_string(ds_regrid.time.values[-1],
                                             unit='D').replace('-', '')

            # Define the output file name based on the start and end dates and
            # the variable name
            var_name = os.path.basename(var_exp_dir)
            output_filename = f"{var_name}_{start_date}_{end_date}.nc"
            output_path = os.path.join(output_dir, output_filename)
            # Print the output file path
            print("Output file:", output_path)

            # Save the merged and regridded dataset to a new netCDF file
            ds_regrid.to_netcdf(output_path)

            # Print the time elapsed for the merging process
            elapsed_time = datetime.now() - now
            print("Elapsed time:", elapsed_time)
        except Exception as e:
            print(f"Error occurred while processing files: {e}")

