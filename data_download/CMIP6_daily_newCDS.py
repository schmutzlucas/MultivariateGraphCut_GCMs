from concurrent.futures import ThreadPoolExecutor
import cdsapi
import os
import time

# Initialize the CDS API client
client = cdsapi.Client(wait_until_complete=True, delete=True)

# Define model list from a file
with open('model_list.txt', 'r') as file:
    model_list = file.read().splitlines()

# Define variables and their short names for CMIP6
variable_dict = {
    'near_surface_air_temperature': 'tas',
    'precipitation': 'pr',
    'sea_level_pressure': 'psl',
    'near_surface_wind_speed': 'sfcWind'
}

# Define the experiment and year range
experiment = 'historical'
start_year = 1850
end_year = 1949

# Define the full date range for the filename
date_range = f"{start_year}0101-{end_year}1231"

# Function to download all data in one request for a given model and variable
def cds_api_call(model, variable, shortname, save_dir):
    # Construct the target filename and a temporary download filename
    target = f"{save_dir}/{shortname}_CMIP6_{experiment}_{model}_{date_range}.zip"
    temp_target = f"{target}.part"

    # Check if the final target file already exists
    if os.path.exists(target):
        print(f"File already exists: {target}")
        return

    # If a temporary file exists, it indicates an incomplete download
    if os.path.exists(temp_target):
        print(f"Resuming download for {temp_target}...")
    else:
        print(f"Starting download: {target}")

    # Define the request parameters
    request = {
        'temporal_resolution': 'daily',
        'experiment': experiment,
        'variable': variable,
        'model': model,
        'year': [str(y) for y in range(start_year, end_year + 1)],
        'month': [f"{m:02d}" for m in range(1, 13)],  # All months
        'day': [f"{d:02d}" for d in range(1, 32)]  # All days
    }

    # Attempt download with retries and failsafe mechanism
    max_retries = 3
    for attempt in range(max_retries):
        try:
            client.retrieve('projections-cmip6', request, temp_target)
            os.rename(temp_target, target)  # Rename to final target upon success
            print(f"Downloaded: {target}")
            return  # Exit if download is successful

        except Exception as e:
            print(f"Attempt {attempt + 1} failed for {target}: {e}")
            time.sleep(2)  # Short delay before retry

    print(f"Failed to download {target} after {max_retries} attempts.")

# Define the main function to run the data retrieval
def main():
    # Base directory to save CMIP6 files
    base_save_dir = os.path.join(os.getcwd(), "CMIP6_data")

    # Create a thread pool for concurrent downloads
    with ThreadPoolExecutor(max_workers=16) as executor:
        for model in model_list:
            for variable, shortname in variable_dict.items():
                # Create directory for each model and variable
                save_dir = f"{base_save_dir}/{experiment}/{model}/{shortname}"
                os.makedirs(save_dir, exist_ok=True)

                # Submit download task for the entire date range
                executor.submit(cds_api_call, model, variable, shortname, save_dir)
                # Small wait time between requests to avoid overloading the server
                time.sleep(0.001)

    print("Download process completed!")

# Execute main function
if __name__ == "__main__":
    main()