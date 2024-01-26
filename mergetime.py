import os
import shutil
import tempfile
import xarray as xr

# Define the root directory containing all the model folders
root_dir = 'Y:\LucasSchmutz\MultivariateGraphCut_GCMs\download_day_unzip'

# Define the directory where the merged files will be stored
merged_dir = 'Y:\LucasSchmutz\MultivariateGraphCut_GCMs\download_day_merged'

# Define the temporary directory to store the merged files before writing to network
temp_dir = tempfile.mkdtemp()

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
                temp_file_path = os.path.join(temp_dir, filename)
                shutil.copy(input_file_path, temp_file_path)
                input_files.append(temp_file_path)

        # Print the input files being merged
        print(f"Merging {len(input_files)} files:")
        print('\n'.join(input_files))

        # Merge the netCDF files using xarray
        ds = xr.open_mfdataset(input_files, combine='nested', concat_dim='time')

        # Define output file name for merged file
        output_filename = os.path.basename(input_files[0]).replace("_19500101-",
                                                                   "_").replace(
            "_v", "_merged_v")

        # Write merged file to local temp directory
        temp_merged_path = os.path.join(temp_dir, output_filename)
        ds.to_netcdf(temp_merged_path)

        # Write merged file to network folder
        output_path = os.path.join(output_dir, output_filename)
        shutil.copy(temp_merged_path, output_path)

        # Remove the temporary input files
        for temp_file_path in input_files:
            os.remove(temp_file_path)

        # Remove the temporary merged file
        os.remove(temp_merged_path)

# Remove the temporary directory
os.rmdir(temp_dir)
