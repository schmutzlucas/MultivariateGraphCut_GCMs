import os
import zipfile
import time
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(
    filename='bad_zip_log.txt',
    level=logging.INFO,
    format='%(asctime)s:%(levelname)s:%(message)s'
)

def unzip_folders(base_path, output_base_path):
    """
    Unzips all zip files within folders and subfolders starting from the base path.
    """
    for root, _, files in os.walk(base_path):  # Walk through all directories and files
        relative_path = os.path.relpath(root, base_path)
        output_path = os.path.join(output_base_path, relative_path)
        os.makedirs(output_path, exist_ok=True)

        for filename in files:
            if filename.endswith('.zip'):
                filepath = os.path.join(root, filename)
                unzip_file(filepath, output_path)

def unzip_file(filepath, output_path):
    """
    Unzips a single zip file into the specified output path.
    """
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{current_time}] Starting to unzip: {filepath}")
    logging.info(f"Starting to unzip: {filepath}")
    start_time = time.time()

    try:
        with zipfile.ZipFile(filepath, 'r') as zip_ref:
            zip_ref.extractall(output_path)  # Extract files to the target output directory
    except zipfile.BadZipFile:
        print(f"Bad zip file: {filepath}")
        logging.error(f"Bad zip file: {filepath}")
    except Exception as e:
        print(f"Error unzipping {filepath}: {e}")
        logging.error(f"Error unzipping {filepath}: {e}")
    else:
        end_time = time.time()
        print(f"Finished unzipping: {filepath}")
        print(f"Time taken: {end_time - start_time:.2f} seconds")
        logging.info(f"Finished unzipping: {filepath}")
        logging.info(f"Time taken: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    # Define paths
    path_to_folders = 'data_download/CMIP6_data'
    output_base_path = 'download_day_unzip_new'

    # Start unzipping process
    unzip_folders(path_to_folders, output_base_path)
