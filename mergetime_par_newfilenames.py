import os
import xarray as xr
import numpy as np
from datetime import datetime
import functools
import multiprocessing
import dask
import netCDF4
import h5netcdf
import scipy


def process_model_dir(root_dir, model_dir, merged_dir):
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
        ds = xr.open_mfdataset(input_files, combine='nested', concat_dim='time')

        grid_path = 'my_grid.txt'
        with open(grid_path, 'r') as f:
            grid_info = f.read()
        grid_dict = {}
        for line in grid_info.split('\n'):
            if '=' in line:
                key, val = line.split('=')
                grid_dict[key.strip()] = float(val)

        # Define the output grid
        output_lon = np.linspace(grid_dict['xfirst'],
                                 grid_dict['xfirst'] + grid_dict['xinc'] * (
                                         grid_dict['xsize'] - 1),
                                 grid_dict['xsize'])
        output_lat = np.linspace(grid_dict['yfirst'],
                                 grid_dict['yfirst'] + grid_dict['yinc'] * (
                                         grid_dict['ysize'] - 1),
                                 grid_dict['ysize'])

        # Define the interpolation method (e.g. 'linear', 'nearest', etc.)
        interp_method = 'linear'

        # Merge and regrid the netCDF files using xarray
        ds = xr.open_mfdataset(input_files, combine='nested', concat_dim='time')
        ds_regrid = ds.interp(lon=output_lon, lat=output_lat,
                              method=interp_method)

        # Extract the start and end dates from the input files and convert
        # them to datetime objects
        start_date = ds_regrid.time.values[0].strftime('%Y%m%d')
        end_date = ds_regrid.time.values[-1].strftime('%Y%m%d')

        # Define output file name for merged file
        input_filename = os.path.basename(input_files[0])
        var = input_filename.split("_")[0]
        model = input_filename.split("_")[2]
        exp = input_filename.split("_")[3]
        ens = input_filename.split("_")[4]
        version = input_filename.split("_")[-1].split(".")[0]
        output_filename = f"{var}_{model}_{exp}_{ens}_{start_date}-{end_date}_{version}.nc"

        # Write merged and regridded file to network folder
        output_path = os.path.join(output_dir, output_filename)
        ds_regrid.to_netcdf(output_path)


if __name__ == '__main__':
    # Define the root directory containing all the model folders
    root_dir = 'Y:\LucasSchmutz\MultivariateGraphCut_GCMs\download_day_unzip'

    # Define the directory where the merged files will be stored
    merged_dir = 'Y:\LucasSchmutz\MultivariateGraphCut_GCMs\download_day_merged'

    # Define the number of processes to use
    num_processes = 64

    # Create a list of model directories
    model_dirs = [d for d in os.listdir(root_dir) if
                  os.path.isdir(os.path.join(root_dir, d))]

    # Create a partial function to pass the fixed arguments to
    # process_model_dir()
    partial_func = functools.partial(process_model_dir, root_dir,
                                     merged_dir=merged_dir)

    # Use multiprocessing.Pool to parallelize the execution of
    # process_model_dir()
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(partial_func, model_dirs)
