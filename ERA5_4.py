import cdsapi
import os

# Create a CDS API client object
c = cdsapi.Client()

# Define the list of variables to retrieve
variable_list = ['2m_temperature']

# Define the list of years to retrieve data for
YEARS = ["1968", "1969", "1970"]

# Define the list of months to retrieve data for
MONTHS = [
    "01", "02", "03", "04", "05", "06",
    "07", "08", "09", "10", "11", "12"
]

# Loop over each year, month, and variable
for year in YEARS:
    for month in MONTHS:
        for variable in variable_list:
            # Set the output filename for the downloaded file
            name = f"ERA5/tas_daily_{year}_{month}.nc"
            # Check if the file already exists before downloading it
            if not os.path.exists(name):
                # Call the CDS API and retrieve data
                result = c.service(
                    "tool.toolbox.orchestrator.run_workflow.user-apps.app-c3s-daily-era5-statistics",
                    params={
                        "code": "",
                        "kwargs": {
                            "area": {
                                "lat": [-90, 90],
                                "lon": [-180, 180]
                            },
                            "dataset": "reanalysis-era5-single-levels",
                            "frequency": "1-hourly",
                            "grid_e5": "1.0/1.0",
                            "month": month,
                            "pressure_level_e5sl": "-",
                            "product_type": "reanalysis",
                            "statistic": "daily_mean",
                            "time_zone": "UTC+00:00",
                            "variable_e5sl": variable,
                            "year_e5sl": year
                        },
                        "workflow_name": "application"
                    },
                )
                c.download(result)

print('finished')
