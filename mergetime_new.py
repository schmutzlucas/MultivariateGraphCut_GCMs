import os
import shutil
import xarray as xr

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
        start_dates = [os.path.basename(input_file).split('_')[3] for input_file in input_files]
        end_dates = [os.path.basename(input_file).split('_')[4][:8] for input_file in input_files]
        start_date = min(start_dates)
        end_date = max(end_dates)

        # Define output file name for merged file
        output_filename = os.path.basename(input_files[0]).replace(start_date + "-",
                                                                   "").replace(
            "_v", f"_{start_date}-{end_date}_v")

        # Write merged file to network folder
        output_path = os.path.join(output_dir, output_filename)
        ds.to_netcdf(output_path)

# Remove the temporary directory if it exists
if os.path.exists(temp_dir):
    os.rmdir(temp_dir)
