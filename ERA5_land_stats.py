from concurrent.futures import ThreadPoolExecutor
import cdsapi
import os
import time


# Create a CDS API client object
c = cdsapi.Client()

# Define a function to call the CDS API and retrieve data
def cds_api_call(year, month, variable, name):
    result = c.service(
        "tool.toolbox.orchestrator.workflow",
        params={
            "realm": "user-apps",
            "project": "app-c3s-daily-era5-statistics",
            "version": "master",
            "kwargs": {
                "area": {
                    "lat": [-90, 90],
                    "lon": [-180, 180]
                },
                "dataset": "reanalysis-era5-land",
                "frequency": "1-hourly",
                "grid_e5": "0.1/0.1",
                "month": month,
                "pressure_level_e5sl": "-",
                "product_type": "reanalysis",
                "statistic": "daily_mean",
                "time_zone": "UTC+00:00",
                "variable_e5sl": variable,
                "year_e5sl": year
            },
            "workflow_name": "application"
        })


    # Find the filename in the result
    temp_filename = None
    for item in result:
        if 'location' in item: \
                temp_filename = item['location'].rsplit('/', 1)[-1]
        break
    print('\n')
    print(temp_filename)
    print('\n')
    if temp_filename:
        # Download the file with the temporary filename
        c.download(result)

        # After downloading, rename it
        os.rename(temp_filename, name)
    else:
        print("Could not find filename in API response")


# Define the list of variables to retrieve
# Map ERA5 variables to CMIP6
variable_map = {
    '2m_temperature': 'tas',
    'total_precipitation': 'pr',
    'surface_pressure' : 'sp'
    # 'mean_sea_level_pressure': 'psl'
}

# Define the list of years to retrieve data for
# Define the list of years to retrieve data for
YEARS = list(map(str, range(1958, 1970)))


# Define the list of months to retrieve data for
MONTHS = [
    "01", "02", "03", "04", "05", "06",
    "07", "08", "09", "10", "11", "12"
]

# Define the main function to run the data retrieval
def main():
    # Ensure the necessary directory exists
    os.makedirs('data/ERA5_land/', exist_ok=True)
    # create a thread pool with n worker threads
    with ThreadPoolExecutor(max_workers=10) as exe:
        for year in YEARS:
            for month in MONTHS:
                # Retrieve data for each variable
                for era5_var, cmip6_var in variable_map.items():
                    # Wait 1 second before the next API call
                    # Needed because faster calls timeout the threads
                    time.sleep(1)
                    # Set the output filename for the downloaded file
                    os.makedirs(f'data/ERA5_land/{cmip6_var}/', exist_ok=True)
                    name = f"data/ERA5_land/{cmip6_var}/{cmip6_var}_ERA5_{year}{month}01-{year}{month}31.nc"
                    # Check if the file already exists before downloading it
                    if not os.path.exists(name):
                        try:
                            # Submit a job to the thread pool to download the
                            # data
                            exe.submit(cds_api_call, year, month, era5_var, name)
                        finally:
                            pass
    print('finished')

# Call the main function to start the data retrieval
main()



