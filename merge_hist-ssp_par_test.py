import logging
import os
import subprocess
import glob
import threading
from concurrent.futures import ThreadPoolExecutor

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(threadName)s - %(message)s')

def extract_dates(filename):
    parts = filename.split("_")[::-1]
    for part in parts:
        if "-" in part and part.split("-")[0].isdigit() and len(part.split("-")[0]) == 8:
            start_date = part.split("-")[0][:8]
            end_date = part.split("-")[1][:8]
            return start_date, end_date
    return None, None

def process_group(model_dir, var_dir, input_files, output_base_dir):
    thread_name = threading.current_thread().name
    logging.info(f"{thread_name}: Starting task for {model_dir}/{var_dir} with {len(input_files)} files")

    if len(input_files) > 1:
        dates = [extract_dates(file) for file in input_files]
        start_date = min(date[0] for date in dates if date[0] is not None)
        end_date = max(date[1] for date in dates if date[1] is not None)

        input_filename = os.path.basename(input_files[0])
        var = input_filename.split("_")[0]
        model = input_filename.split("_")[1]
        output_filename = f"{var}_{model}_{start_date}-{end_date}.nc"
        output_dir = os.path.join(output_base_dir, model_dir, var_dir)
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, output_filename)

        cdo_command = ["cdo", "-O", "mergetime"] + input_files + [output_path]
        subprocess.run(cdo_command)
        logging.info(f"{thread_name}: Finished task for {model_dir}/{var_dir}")
    else:
        logging.info(f"{thread_name}: Skipping merging for {model_dir}/{var_dir}. Only one file found.")


root_dir = '/mnt/y/LucasSchmutz/MultivariateGraphCut_GCMs/merged_regridded_par_new/'
output_base_dir = '/mnt/y/LucasSchmutz/MultivariateGraphCut_GCMs/data/CMIP6_merged_all_new/'

tasks = []

for model_dir in os.listdir(root_dir):
    model_path = os.path.join(root_dir, model_dir)
    for var_dir in os.listdir(model_path):
        var_path = os.path.join(model_path, var_dir)
        nc_files = glob.glob(os.path.join(var_path, "*.nc"))

        if nc_files:  # Proceed only if there are files
            grouped_files = {}
            for file in nc_files:
                parts = os.path.basename(file).split("_")
                key = "_".join([parts[0], parts[1]])  # variable, model
                if key not in grouped_files:
                    grouped_files[key] = []
                grouped_files[key].append(file)

            for key, input_files in grouped_files.items():
                tasks.append((model_dir, var_dir, input_files, output_base_dir))

num_workers = 22  # Adjust based on your system capabilities

with ThreadPoolExecutor(max_workers=num_workers) as executor:
    futures = [executor.submit(process_group, *task) for task in tasks]

    for future in futures:
        future.result()  # Wait for all tasks to complete
