import os

base_dir = '/mnt/y/LucasSchmutz/MultivariateGraphCut_GCMs/data/CMIP6_merged/'
new_end_date = "21001231"

# Iterate through all directories and files
for root, dirs, files in os.walk(base_dir):
    for file in files:
        # If it's a .nc file
        if file.endswith('.nc'):
            # Extract information from the path
            path_parts = os.path.normpath(root).split(os.sep)
            model = path_parts[-2]
            var = path_parts[-1]

            # Extract start date from the filename
            date_range = file.split("_")[2]
            start_date = date_range.split('-')[0]

            # Define new filename with new end date
            new_filename = f"{var}_{model}_{start_date}-{new_end_date}_merged.nc"

            # Define full old and new paths
            old_path = os.path.join(root, file)
            new_path = os.path.join(root, new_filename)

            # Rename the file
            os.rename(old_path, new_path)
