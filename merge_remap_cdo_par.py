import os
import subprocess
import glob
import datetime
import logging
from concurrent.futures import ThreadPoolExecutor

root_dir = 'Y:\\LucasSchmutz\\MultivariateGraphCut_GCMs\\download_day_unzip'
output_dir = 'Y:\\LucasSchmutz\\MultivariateGraphCut_GCMs\\merged_regridded'
grid_file = 'Y:\\LucasSchmutz\\MultivariateGraphCut_GCMs\\my_grid.txt'

logging.basicConfig(filename='merge_remap_cdo.log', level=logging.INFO)

def process_files(input_files, output_path):
    try:
        print(f"Merging {len(input_files)} files:")
        for file in input_files:
            print(file)

        cdo_command = ["cdo", "-O", "remapbil," + grid_file,
                       "-mergetime"] + input_files + [output_path]

        subprocess.run(cdo_command)
    except Exception as e:
        logging.error(f'Error while processing {input_files}: {e}')

tasks = []

for model_dir in os.listdir(root_dir):
    model_path = os.path.join(root_dir, model_dir)

    for var_dir in os.listdir(model_path):
        var_path = os.path.join(model_path, var_dir)

        nc_files = glob.glob(os.path.join(var_path, "*.nc"))

        grouped_files = {}
        for file in nc_files:
            key = "_".join(os.path.basename(file).split("_")[:-2])
            if key not in grouped_files:
                grouped_files[key] = []
            grouped_files[key].append(file)

        for key, input_files in grouped_files.items():
            start_date = input_files[0].split("_")[-2].split("-")[0]
            end_date = input_files[-1].split("_")[-2].split("-")[1]

            input_filename = os.path.basename(input_files[0])
            var = input_filename.split("_")[0]
            model = input_filename.split("_")[2]
            exp = input_filename.split("_")[3]
            ens = input_filename.split("_")[4]
            version = input_filename.split("_")[-1].split(".")[0]
            output_filename = f"{var}_{model}_{exp}_{ens}_{start_date}-{end_date}_merged_regridded_{version}.nc"
            output_path = os.path.join(output_dir, output_filename)

            tasks.append((input_files, output_path))

max_workers = 4  # Change this value to set the number of parallel tasks
with ThreadPoolExecutor(max_workers=max_workers) as executor:
    for task in tasks:
        executor.submit(process_files, task[0], task[1])
