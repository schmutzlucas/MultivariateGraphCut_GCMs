import os
import xarray as xr
import pandas as pd
from datetime import datetime


def merge_files(base_dir):
    print(f"Current time: {datetime.now()}")

    # Discover variable folders
    variable_folders = [os.path.join(base_dir, d) for d in os.listdir(base_dir)
                        if os.path.isdir(os.path.join(base_dir, d))]

    for variable_folder in variable_folders:
        file_list = []
        for root, _, files in os.walk(variable_folder):
            for file in files:
                if file.endswith(".nc"):
                    file_list.append(os.path.join(root, file))

        if file_list:
            try:
                print(f"Merging {len(file_list)} files:")
                for file in file_list:
                    print(file)
                merged_data = xr.combine_by_coords(
                    [xr.open_dataset(f) for f in file_list],
                    combine_attrs='override')
                output_file = os.path.join(variable_folder,
                                           f"merged_{os.path.basename(variable_folder)}.nc")
                print(f"Saving merged data to {output_file}")
                merged_data.to_netcdf(output_file)
                print(f"Saved {output_file}")
            except Exception as e:
                print(f"An error occurred while processing {file_list}: {e}")


if __name__ == "__main__":
    base_dir = "Y:/LucasSchmutz/MultivariateGraphCut_GCMs/download_day_unzip"
    merge_files(base_dir)
