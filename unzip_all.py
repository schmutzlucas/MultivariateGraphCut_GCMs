import os
import zipfile
import time
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(filename='bad_zip_log.txt', level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')

def unzip_folders(path):
    for foldername in os.listdir(path):
        folderpath = os.path.join(path, foldername)
        if os.path.isdir(folderpath):
            unzip_folder(folderpath)

def unzip_folder(path):
    output_path = os.path.join('download_day_unzip_new', os.path.basename(path))
    os.makedirs(output_path, exist_ok=True)
    for filename in os.listdir(path):
        filepath = os.path.join(path, filename)
        if os.path.isfile(filepath) and filename.endswith('.zip'):
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[{current_time}] Starting to unzip: {filepath}")
            logging.info(f"Starting to unzip: {filepath}")
            start_time = time.time()
            try:
                with zipfile.ZipFile(filepath, 'r') as zip_ref:
                    zip_ref.extractall(os.path.join(output_path, filename[:-4]))
            except zipfile.BadZipFile:
                print(f"Bad zip file: {filepath}")
                logging.error(f"Bad zip file: {filepath}")
            else:
                end_time = time.time()
                print(f"Finished unzipping: {filepath}")
                print(f"Time taken: {end_time - start_time:.2f} seconds")
                logging.info(f"Finished unzipping: {filepath}")
                logging.info(f"Time taken: {end_time - start_time:.2f} seconds")

path_to_folders = 'download_day'
unzip_folders(path_to_folders)
