import os
import shutil
import xarray as xr
import numpy as np
from datetime import datetime
from scipy.interpolate import griddata

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

# Define the destination grid for regridding
grid_x, grid_y = np.meshgrid(np.linspace(xfirst, xfirst + xinc * (xsize - 1), xsize),
                             np.linspace(yfirst, yfirst + yinc * (ysize - 1), ysize))

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
            new_data = []
            for var in ds.data_vars:
                data = ds[var].values
                old_lat = ds[var].lat.values
                old_lon = ds[var].lon.values
                points = np.array([old_lon.ravel(), old_lat.ravel()]).T
                new_points = np.array([grid_x.ravel(), grid_y.ravel()]).T
                new_data.append(griddata(points, data.ravel(), new_points,
                                         method='linear').reshape(ysize, xsize))

            # Create a new xarray dataset with the regridded data
            new_ds = xr.Dataset(
                {var: (['time', 'lat', 'lon'], new_data[i]) for i, var in
                 enumerate(ds.data_vars)})
            new_ds['time'] = ds['time']
            new_ds['lat'] = (['lat'], grid_y[:, 0])
            new_ds['lon'] = (['lon'], grid_x[0, :])

            # Write the new dataset to a netCDF file
            output_file_path = os.path.join(output_dir, f"{var_exp_dir}.nc")
            new_ds.to_netcdf(output_file_path)

            # Print a success message
            print(f"Regridded file saved to {output_file_path}")
        except Exception as e:
            print(f"Error occurred while regridding: {e}")
            continue
