import os
import subprocess
import glob
import datetime

root_dir = '/home/lschmutz1/LucasSchmutz/MultivariateGraphCut_GCMs/download_day_unzip/'
output_base_dir = '/home/lschmutz1/LucasSchmutz/MultivariateGraphCut_GCMs/merged_regridded/'
grid_file = '/home/lschmutz1/LucasSchmutz/MultivariateGraphCut_GCMs/my_grid.txt'

# Iterate through each model folder
for model_dir in os.listdir(root_dir):
    model_path = os.path.join(root_dir, model_dir)

    # Iterate through each variable folder
    for var_dir in os.listdir(model_path):
        var_path = os.path.join(model_path, var_dir)

        # Find all NetCDF files in the variable folder
        nc_files = glob.glob(os.path.join(var_path, "*.nc"))

        # Group input files by variable, model, experiment and ensemble member
        grouped_files = {}
        for file in nc_files:
            key = "_".join(os.path.basename(file).split("_")[:-2])
            if key not in grouped_files:
                grouped_files[key] = []
            grouped_files[key].append(file)

        # Process each group of input files
        for key, input_files in grouped_files.items():
            print(f"Merging {len(input_files)} files:")
            for file in input_files:
                print(file)

            # Get the start and end dates from the input files
            start_date = input_files[0].split("_")[-2].split("-")[0]
            end_date = input_files[-1].split("_")[-2].split("-")[1]

            # Define output file name for merged and regridded file
            input_filename = os.path.basename(input_files[0])
            var = input_filename.split("_")[0]
            model = input_filename.split("_")[2]
            exp = input_filename.split("_")[3]
            ens = input_filename.split("_")[4]
            version = input_filename.split("_")[-1].split(".")[0]
            output_filename = f"{var}_{model}_{exp}_{ens}_{start_date}-{end_date}_merged_regridded_{version}.nc"

            # Create output directory with model_dir and var_dir structure
            output_dir = os.path.join(output_base_dir, model_dir, var_dir)
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, output_filename)

            # Define CDO command
            cdo_command = ["cdo", "-O", "remapbil," + grid_file,
                           "-mergetime"] + input_files + [output_path]

            # Execute the CDO command using subprocess
            subprocess.run(cdo_command)
