import os
import subprocess
from concurrent.futures import ThreadPoolExecutor


def process_file(file_path, output_file_path):
    # CDO command to multiply data
    cdo_command = f"cdo mulc,86400 {file_path} {output_file_path}"
    try:
        # Execute the CDO command
        subprocess.run(cdo_command, check=True, shell=True)
        print(f"Processed {file_path} -> {output_file_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Failed to process {file_path}: {str(e)}")
        return False


def process_files(root_dir, output_root):
    # Create the output root directory if it doesn't exist
    os.makedirs(output_root, exist_ok=True)

    # Prepare tasks for the ThreadPoolExecutor
    tasks = []
    for subdir, dirs, files in os.walk(root_dir):
        for file in files:
            if file.startswith("pr_") and file.endswith(".nc"):
                file_path = os.path.join(subdir, file)
                new_subdir = subdir.replace(root_dir, output_root)
                os.makedirs(new_subdir, exist_ok=True)
                output_file_path = os.path.join(new_subdir, file)
                tasks.append((file_path, output_file_path))

    # Number of threads to use; you can adjust this number as needed
    num_threads = 4

    # Process files in parallel using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = executor.map(lambda x: process_file(*x), tasks)
        for result in results:
            if not result:
                print("Some files failed to process.")


# Example usage
root_dir = 'data/CMIP6_merged_all'
output_root = 'data/CMIP6_corrected'
process_files(root_dir, output_root)
