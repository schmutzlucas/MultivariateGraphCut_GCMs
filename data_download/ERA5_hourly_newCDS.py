from concurrent.futures import ThreadPoolExecutor
import cdsapi
import os
import time

# Create a CDS API client object
client = cdsapi.Client()

# Define the list of variables with long and short names for ERA5
variable_dict = {
    '2m_temperature': 'tas',
    'total_precipitation': 'pr',
    'mean_sea_level_pressure': 'msl'
}

# Define a function to delete completed requests
def delete_completed_request(request_id):
    try:
        client.delete(request_id)
        print(f"Deleted request: {request_id}")
    except Exception as e:
        print(f"Failed to delete request {request_id}: {e}")

# Define a function to call the CDS API and retrieve data
def cds_api_call(year, month, variable, shortname, save_dir):
    # Construct the target filename and a temporary download filename
    target = f"{save_dir}/{shortname}_ERA5_{year}{month:02d}.nc"
    temp_target = f"{target}.part"

    # Define the request parameters
    request = {
        "product_type": "reanalysis",
        "variable": [variable],
        "year": str(year),
        "month": [f"{month:02d}"],
        "day": [f"{d:02d}" for d in range(1, 32)],  # All days in the month
        "time": [  # All hours of the day
            "00:00", "01:00", "02:00", "03:00", "04:00", "05:00",
            "06:00", "07:00", "08:00", "09:00", "10:00", "11:00",
            "12:00", "13:00", "14:00", "15:00", "16:00", "17:00",
            "18:00", "19:00", "20:00", "21:00", "22:00", "23:00"
        ],
        "format": "netcdf",  # Output format specified as NetCDF
        "target": "disk"  # Retrieve from disk if available
    }

    # Check if the final target file already exists
    if os.path.exists(target):
        print(f"File already exists: {target}")
        return

    # If a temporary file exists, it indicates an incomplete download
    if os.path.exists(temp_target):
        print(f"Resuming download for {temp_target}...")
    else:
        print(f"Starting download: {target}")

    # Download the file to the temporary target first
    try:
        result = client.retrieve('reanalysis-era5-single-levels', request, temp_target)
        # Rename to the final target file name upon successful download
        os.rename(temp_target, target)
        print(f"Downloaded: {target}")

        # Delete the completed request from the server
        request_id = result.get("request_id", None)
        if request_id:
            delete_completed_request(request_id)

    except Exception as e:
        print(f"Error downloading {target}: {e}")
        # If there's an error, the partially downloaded file is left as is

# Define the main function to run the data retrieval
def main():
    # Set the year range for data retrieval
    start_year = 1962
    end_year = 1970

    # Define the base directory for saving files
    base_save_dir = os.path.join(os.getcwd(), "ERA5")

    # Create a thread pool with a defined number of worker threads
    with ThreadPoolExecutor(max_workers=32) as executor:
        # Iterate through each variable and its short name
        for longname, shortname in variable_dict.items():
            # Create a save directory for each variable if it doesn't exist
            save_dir = f"{base_save_dir}/{shortname}"
            os.makedirs(save_dir, exist_ok=True)

            # Iterate through years and months to define the download period
            for year in range(start_year, end_year + 1):
                for month in range(1, 13):  # Loop through all months (1-12)
                    # Submit a job to the thread pool to download data for each month
                    executor.submit(cds_api_call, year, month, longname, shortname, save_dir)
                    # Add a small wait time between requests to avoid overloading the server
                    time.sleep(0.3)

    print('Download completed!')

# Call the main function to start the data retrieval
if __name__ == "__main__":
    main()
