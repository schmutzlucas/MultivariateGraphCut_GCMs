import xarray as xr
import os
import glob

root_dir = '/mnt/y/LucasSchmutz/MultivariateGraphCut_GCMs/data/CMIP6/'
output_base_dir = '/mnt/y/LucasSchmutz/MultivariateGraphCut_GCMs/data/CMIP6_merged/'

# Helper function to extract start and end dates from the combined array
def extract_dates(combined):
    start_date = combined['time'].min().dt.strftime('%Y%m%d').item()
    end_date = combined['time'].max().dt.strftime('%Y%m%d').item()
    return start_date, end_date

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
            key = "_".join(parts[:2])  # variable, model
            if key not in grouped_files:
                grouped_files[key] = []
            grouped_files[key].append(file)

        print(f"Grouped files: {grouped_files}")

        # Process each group of input files
        for key, input_files in grouped_files.items():
            print(f"Merging {len(input_files)} files:")
            for file in input_files:
                print(file)

            # Merge files
            combined = None
            datasets = []  # Store the opened datasets
            try:
                for file in input_files:
                    print(f"Processing file: {file}")

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

                # Get the start and end dates from the combined array
                start_date, end_date = extract_dates(combined)

                # Define output file name for merged file
                input_filename = os.path.basename(input_files[0])
                var = input_filename.split("_")[0]
                model = input_filename.split("_")[1]
                output_filename = f"{var}_{model}_{start_date}-{end_date}_merged.nc"

                # Create output directory with model_dir and var_dir structure
                output_dir = os.path.join(output_base_dir, model_dir, var_dir)
                os.makedirs(output_dir, exist_ok=True)
                output_path = os.path.join(output_dir, output_filename)

                # Save the combined dataset
                print(f"Saving merged file: {output_filename}")
                combined.to_netcdf(output_path)

                print(f"Merge complete: {output_filename}")

            finally:
                # Close the datasets
                for dataset in datasets:
                    dataset.close()


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
            print(f"Saving merged file: {output_filename}")
            combined.to_netcdf(output_path)
            combined.close()  # Close the combined dataset

            print(f"Merge complete: {output_filename}")
