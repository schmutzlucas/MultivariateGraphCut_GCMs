import xarray as xr
import os
import glob
import dask

root_dir = '/mnt/y/LucasSchmutz/MultivariateGraphCut_GCMs/data/CMIP6/'
output_base_dir = '/mnt/y/LucasSchmutz/MultivariateGraphCut_GCMs/data/CMIP6_merged/'


def extract_dates(filename):
    parts = filename.split("_")[::-1]
    for part in parts:
        if "-" in part and part.split("-")[0].isdigit() and len(part.split("-")[0]) == 8:
            return part.split("-")[0][:4], part.split("-")[1][:4]
    return None, None


def process_files(input_files):
    print(f"Merging {len(input_files)} files:")
    for file in input_files:
        print(file)

    datasets = [dask.delayed(xr.open_dataset)(f, chunks={}) for f in input_files]
    combined = xr.concat(dask.compute(*datasets), dim='time')
    combined = combined.sortby('time')

    dates = [extract_dates(file) for file in input_files]
    start_date = min(date[0] for date in dates if date[0] is not None)
    end_date = max(date[1] for date in dates if date[1] is not None)

    input_filename = os.path.basename(input_files[0])
    var = input_filename.split("_")[0]
    model = input_filename.split("_")[2]
    ens = input_filename.split("_")[4]
    version = input_filename.split("_")[-1].split(".")[0]
    output_filename = f"{var}_{model}_{ens}_{start_date}-{end_date}_merged_{version}.nc"

    output_dir = os.path.join(output_base_dir, model_dir, var_dir)
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_filename)

    combined.to_netcdf(output_path, engine='h5netcdf')

    for ds in datasets:
        ds[0].close()  # Need to access the computed xarray Dataset, hence the [0]
    combined.close()


for model_dir in os.listdir(root_dir):
    model_path = os.path.join(root_dir, model_dir)
    for var_dir in os.listdir(model_path):
        var_path = os.path.join(model_path, var_dir)
        nc_files = glob.glob(os.path.join(var_path, "*.nc"))
        print(f"Found {len(nc_files)} .nc files: {nc_files}")

        grouped_files = {}
        for file in nc_files:
            parts = os.path.basename(file).split("_")
            key = "_".join([parts[0], parts[1]])
            if key not in grouped_files:
                grouped_files[key] = []
            grouped_files[key].append(file)

        print(f"Grouped files: {grouped_files}")

        for key, input_files in grouped_files.items():
            process_files(input_files)
