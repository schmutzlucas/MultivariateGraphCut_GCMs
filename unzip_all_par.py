import os
import zipfile
from concurrent.futures import ThreadPoolExecutor


def unzip_folders(path):
    with ThreadPoolExecutor(max_workers=8) as executor:
        for foldername in os.listdir(path):
            folderpath = os.path.join(path, foldername)
            if os.path.isdir(folderpath):
                executor.submit(unzip_folder, folderpath)


def unzip_folder(path):
    output_path = os.path.join('Y:/LucasSchmutz/MultivariateGraphCut_GCMs'
                               '/download_day_unzip', os.path.basename(path))
    os.makedirs(output_path, exist_ok=True)
    for filename in os.listdir(path):
        filepath = os.path.join(path, filename)
        if os.path.isfile(filepath) and filename.endswith('.zip'):
            with zipfile.ZipFile(filepath, 'r') as zip_ref:
                zip_ref.extractall(os.path.join(output_path, filename[:-4]))


path_to_folders = 'download_day'
unzip_folders(path_to_folders)
