import os
import subprocess
import glob

root_dir = '/mnt/y/LucasSchmutz/MultivariateGraphCut_GCMs/merged_regridded_par/'
output_base_dir = '/mnt/y/LucasSchmutz/MultivariateGraphCut_GCMs/data/CMIP6_merged_all/'

# Helper function to extract date range from filename
def extract_dates(filename):
    # Split the filename and reverse it
    parts = filename.split("_")[::-1]

    # Find the first part that contains a date
    for part in parts:
        if "-" in part and part.split("-")[0].isdigit() and len(part.split("-")[0]) == 8:
            start_date = part.split("-")[0][:8]
            end_date = part.split("-")[1][:8]
            return start_date, end_date
    return None, None

# Iterate through each model folder
for model_dir in os.listdir(root_dir):
    model_path = os.path.join(root_dir, model_dir)

    # Iterate through each variable folder
    for var_dir in os.listdir(model_path):
        var_path = os.path.join(model_path, var_dir)

        # Find all NetCDF files in the variable folder
        nc_files = glob.glob(os.path.join(var_path, "*.nc"))
        print(f"Found {len(nc_files)} .nc files: {nc_files}")

        # Group input files by variable, model, and ensemble member
        grouped_files = {}
        for file in nc_files:
            parts = os.path.basename(file).split("_")
            key = "_".join([parts[0], parts[1]])  # variable, model
            if key not in grouped_files:
                grouped_files[key] = []
            grouped_files[key].append(file)

        print(f"Grouped files: {grouped_files}")

        # Process each group of input files
        for key, input_files in grouped_files.items():
            print(f"Merging {len(input_files)} files:")
            for file in input_files:
                print(file)

            if len(input_files) > 1:
                # Get the start and end dates from the input files
                dates = [extract_dates(file) for file in input_files]
                start_date = min(date[0] for date in dates if date[0] is not None)
                end_date = max(date[1] for date in dates if date[1] is not None)

                # Define output file name for merged file
                input_filename = os.path.basename(input_files[0])
                var = input_filename.split("_")[0]
                model = input_filename.split("_")[1]
                output_filename = f"{var}_{model}_{start_date}-{end_date}.nc"

                # Create output directory with model_dir and var_dir structure
                output_dir = os.path.join(output_base_dir, model_dir, var_dir)
                os.makedirs(output_dir, exist_ok=True)
                output_path = os.path.join(output_dir, output_filename)

                # Define CDO command
                cdo_command = ["cdo", "-O", "mergetime"] + input_files + [output_path]

                # Execute the CDO command using subprocess asynchronously
                subprocess.run(cdo_command)
            else:
                print("Skipping merging. Only one file found.")
