import os
import zipfile
import time
import logging
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

# Set up logging
logging.basicConfig(
    filename='bad_zip_log.txt',
    level=logging.INFO,
    format='%(asctime)s:%(levelname)s:%(message)s'
)

def find_zip_files(base_path):
    """
    Recursively find all zip files in the directory structure.
    """
    zip_files = []
    for root, _, files in os.walk(base_path):
        for filename in files:
            if filename.endswith('.zip'):
                zip_files.append(os.path.join(root, filename))
    return zip_files

def unzip_file(filepath, output_base_path):
    """
    Unzips a single zip file into the corresponding output path.
    """
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{current_time}] Starting to unzip: {filepath}")
    logging.info(f"Starting to unzip: {filepath}")
    start_time = time.time()

    try:
        # Determine output path based on relative path
        relative_path = os.path.relpath(os.path.dirname(filepath), output_base_path)
        output_path = os.path.join(output_base_path, relative_path)
        os.makedirs(output_path, exist_ok=True)

        # Extract the zip file
        with zipfile.ZipFile(filepath, 'r') as zip_ref:
            zip_ref.extractall(output_path)

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

def unzip_folders_parallel(base_path, output_base_path, max_workers=8):
    """
    Unzips all zip files in parallel using ThreadPoolExecutor.
    """
    # Find all zip files in the directory structure
    zip_files = find_zip_files(base_path)

    # Use ThreadPoolExecutor for parallel unzipping
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_file = {
            executor.submit(unzip_file, filepath, output_base_path): filepath for filepath in zip_files
        }

        for future in as_completed(future_to_file):
            filepath = future_to_file[future]
            try:
                future.result()  # This will re-raise exceptions from unzip_file
            except Exception as e:
                print(f"Error processing {filepath}: {e}")
                logging.error(f"Error processing {filepath}: {e}")

if __name__ == "__main__":
    # Define paths
    path_to_folders = 'data_download/CMIP6_data'
    output_base_path = 'data_download/download_day_unzip_new'

    # Start unzipping process with parallelization
    unzip_folders_parallel(path_to_folders, output_base_path, max_workers=8)
