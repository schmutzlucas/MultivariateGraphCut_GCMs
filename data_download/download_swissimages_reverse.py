import csv
import requests
import os
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

# Define the path to the CSV file
csv_file_path = 'ch.swisstopo.swissimage-dop10-d0CZG8Es.csv'

# Define the directory to save the downloaded images
save_directory = 'swissimages'

# Create the directory if it doesn't exist
os.makedirs(save_directory, exist_ok=True)

# Maximum number of retries for a failed download
max_retries = 3
# Timeout for each request (in seconds)
request_timeout = 30
# Timeout for chunk download (in seconds)
chunk_timeout = 10

# Function to extract the filename from the Content-Disposition header
def get_filename_from_response(response):
    if 'Content-Disposition' in response.headers:
        content_disposition = response.headers['Content-Disposition']
        filename = content_disposition.split('filename=')[-1].strip('\"')
        return filename
    else:
        return None

# Function to download an image from a URL with progress bar, size check, and retry mechanism
def download_image(url):
    for attempt in range(max_retries):
        try:
            response = requests.get(url, stream=True, timeout=request_timeout)
            response.raise_for_status()  # Check if the request was successful

            # Extract filename from the response header
            filename = get_filename_from_response(response)
            if not filename:
                filename = os.path.basename(url)

            save_path = os.path.join(save_directory, filename)

            total_size = int(response.headers.get('content-length', 0))

            # Check if the file already exists and has the correct size
            if os.path.exists(save_path):
                existing_file_size = os.path.getsize(save_path)
                if existing_file_size == total_size:
                    print(f"File {filename} already exists and is complete. Skipping download.")
                    return
                else:
                    print(f"File {filename} is incomplete or corrupted (expected {total_size} bytes, found {existing_file_size} bytes). Redownloading.")
                    os.remove(save_path)  # Delete the incomplete or corrupted file

            # Proceed to download the image chunk-by-chunk
            block_size = 1024  # 1 Kibibyte
            t = tqdm(total=total_size, unit='iB', unit_scale=True, desc=filename, leave=False)

            with open(save_path, 'wb') as file:
                for data in response.iter_content(block_size):
                    if not data:
                        raise requests.exceptions.Timeout("Chunk download timed out")
                    t.update(len(data))
                    file.write(data)
                    # Pause briefly between chunks to avoid stalling
                    time.sleep(0.1)
            t.close()

            if total_size != 0 and t.n != total_size:
                print(f"WARNING: {filename} download incomplete.")
            else:
                print(f"Image downloaded: {filename}")
            return  # If download was successful, exit the function

        except (requests.exceptions.RequestException, requests.exceptions.Timeout) as e:
            print(f"Attempt {attempt + 1} failed for {url}. Error: {e}")
            time.sleep(2)  # Wait a bit before retrying

    # If we exit the loop without downloading successfully
    print(f"Failed to download {url} after {max_retries} attempts.")

# Function to download images concurrently
def download_images_concurrently(urls, max_workers=10):
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_url = {executor.submit(download_image, url): url for url in urls}
        for future in tqdm(as_completed(future_to_url), total=len(urls), desc="Overall Progress"):
            try:
                future.result()
            except Exception as e:
                url = future_to_url[future]
                print(f"Error downloading {url}: {e}")

# Read the CSV file and prepare the URLs
urls = []
with open(csv_file_path, mode='r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        image_url = row[0]  # Assuming the image URLs are in the first column
        urls.append(image_url)

# Reverse the list of URLs to start downloading from the end
urls.reverse()

# Download images using multiple workers
download_images_concurrently(urls, max_workers=10)
