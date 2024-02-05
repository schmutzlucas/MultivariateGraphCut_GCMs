import xarray as xr
import os
import glob
from multiprocessing import Pool

root_dir = '/mnt/y/LucasSchmutz/MultivariateGraphCut_GCMs/merged_regridded_par_new'
output_base_dir = '/mnt/y/LucasSchmutz/MultivariateGraphCut_GCMs/data/CMIP6_merged_new/'

# Helper function to extract date range from filename
def extract_dates(filename):
    # Split the filename and reverse it
    parts = filename.split("_")[::-1]

    # Find the first part that contains a date
    for part in parts:
        if "-" in part and part.split("-")[0].isdigit() and len(part.split("-")[0]) == 8:
            return part.split("-")[0][:4], part.split("-")[1][:4]
    return None, None

# Function to merge a group of input files
def merge_files(input_files):
    # Merge files
    combined = None
    datasets = []  # Store the opened datasets
    try:
        for file in input_files:
            # Load the file into an xarray dataset
            dataset = xr.open_dataset(file)
            datasets.append(dataset)

            if combined is None:
                # If it's the first file, initialize the combined dataset
                combined = dataset
            else:
                # Concatenate the current dataset with the combined dataset
                combined = xr.concat([combined, dataset], dim='time')
                combined = combined.sortby('time')  # Ensure time is in ascending order

        # Get the start and end dates from the input files
        dates = [extract_dates(file) for file in input_files]
        start_date = min(date[0] for date in dates if date[0] is not None)
        end_date = max(date[1] for date in dates if date[1] is not None)

        # Define output file name for merged file
        input_filename = os.path.basename(input_files[0])
        var = input_filename.split("_")[0]
        model = input_filename.split("_")[2]
        ens = input_filename.split("_")[4]
        version = input_filename.split("_")[-1].split(".")[0]
        output_filename = f"{var}_{model}_{ens}_{start_date}-{end_date}_merged.nc"

        # Create output directory with model_dir and var_dir structure
        output_dir = os.path.join(output_base_dir, model_dir, var_dir)
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, output_filename)

        # Save the combined dataset
        combined.to_netcdf(output_path)
        combined.close()  # Close the combined dataset

        print(f"Merge complete: {output_filename}")

    finally:
        # Close the datasets
        for dataset in datasets:
            dataset.close()

# Iterate through each model folder
for model_dir in os.listdir(root_dir):
    model_path = os.path.join(root_dir, model_dir)
    print(f"Processing model directory: {model_dir}")

    # Iterate through each variable folder
    for var_dir in os.listdir(model_path):
        var_path = os.path.join(model_path, var_dir)
        print(f"Processing variable directory: {var_dir}")

        # Find all NetCDF files in the variable folder
        nc_files = glob.glob(os.path.join(var_path, "*.nc"))
        print(f"Found {len(nc_files)} .nc files: {nc_files}")

        # Group input files by variable, model, and ensemble member
        grouped_files = {}
        for file in nc_files:
            parts = os.path.basename(file).split("_")
            key = "_".join([parts[0], parts[1]])  # variable, model, ensemble
            if key not in grouped_files:
                grouped_files[key] = []
            grouped_files[key].append(file)

        print(f"Grouped files: {grouped_files}")

        # Create a pool of worker processes
        pool = Pool()

        # Perform merges in parallel
        pool.map(merge_files, grouped_files.values())

        # Close the pool of worker processes
        pool.close()
        pool.join()
