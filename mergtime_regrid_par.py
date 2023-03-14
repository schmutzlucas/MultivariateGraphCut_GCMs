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
import xesmf as xe

def process_model_dir(root_dir, model_dir, merged_dir, regrid_kwargs=None):
    # Create a dask client for parallel processing
    client = dd.Client()

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

        # Merge and regrid the netCDF files using xarray and dask
        dsets = [xr.open_dataset(f) for f in input_files]
        ds = xr.concat(dsets, dim='time', data_vars='minimal', coords='minimal',
                       compat='override')
        ds = ds.chunk({'time': 'auto'})

        if regrid_kwargs is not None:
            # Regrid the data using xesmf and dask


            # Define the target grid for regridding
            target_grid = xr.Dataset({'lat': (['y', 'x'],
                                              np.arange(regrid_kwargs['yfirst'],
                                                        regrid_kwargs[
                                                            'yfirst'] +
                                                        regrid_kwargs['ysize'] *
                                                        regrid_kwargs['yinc'],
                                                        regrid_kwargs['yinc'])),
                                      'lon': (['y', 'x'],
                                              np.arange(regrid_kwargs['xfirst'],
                                                        regrid_kwargs[
                                                            'xfirst'] +
                                                        regrid_kwargs['xsize'] *
                                                        regrid_kwargs['xinc'],
                                                        regrid_kwargs[
                                                            'xinc']))})

            # Create the regridder object
            regridder = xe.Regridder(ds, target_grid, 'bilinear',
                                     periodic=False, reuse_weights=True)

            # Regrid the data
            ds = regridder(ds)

        # Extract the start and end dates from the input files and convert
        # them to datetime objects
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

        if regrid_kwargs is not None:
            # Add suffix to version for regridded files
            version += f"_regrid_{regrid_kwargs['xsize']}x{regrid_kwargs['ysize']}_{regrid_kwargs['gridtype']}"

        output_filename = f"{var}_{model}_{exp}_{ens}_{start_date}-{end_date}_{version}.nc"

        # Write merged file to network folder
        output_path = os.path.join(output_dir, output_filename)
        ds.to_netcdf(output_path, compute=False)
        # Close the dataset to free up resources
        ds.close()

    # Close the dask client
    client.close()


if __name__ == '__main__':
    # Define the root directory containing all the model folders
    root_dir = "/path/to/models"
    # Define the directory where the merged files will be saved
    merged_dir = "/path/to/merged/files"
    # Define the regridding parameters
    regrid_kwargs = {'gridtype': 'lonlat',
                     'xsize': 360,
                     'ysize': 181,
                     'xfirst': 0,
                     'xinc': 1,
                     'yfirst': -90,
                     'yinc': 1}

    # Loop over all model directories in the root directory
    for model_dir in os.listdir(root_dir):
        # Ignore any non-directory files in the root directory
        if not os.path.isdir(os.path.join(root_dir, model_dir)):
            continue

        # Process the model directory
        process_model_dir(root_dir, model_dir, merged_dir, regrid_kwargs)
