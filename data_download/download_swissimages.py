import csv
import requests
import os

# Define the path to the CSV file
csv_file_path = 'data_download/ch.swisstopo.swissimage-dop10-d0CZG8Es.csv'

# Define the directory to save the downloaded images
save_directory = 'data_download/swissimages'

# Create the directory if it doesn't exist
os.makedirs(save_directory, exist_ok=True)

# Function to download an image from a URL
def download_image(url, save_path):
    try:
        response = requests.get(url)
        response.raise_for_status()  # Check if the request was successful
        with open(save_path, 'wb') as file:
            file.write(response.content)
        print(f"Image downloaded: {save_path}")
    except requests.exceptions.RequestException as e:
        print(f"Failed to download {url}. Error: {e}")

# Read the CSV file and download each image
with open(csv_file_path, mode='r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        image_url = row[0]  # Assuming the image URLs are in the first column
        image_name = os.path.basename(image_url)
        save_path = os.path.join(save_directory, image_name)
        download_image(image_url, save_path)
