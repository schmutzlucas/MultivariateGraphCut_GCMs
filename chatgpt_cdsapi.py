import cdsapi
import asyncio


# Define the CDS API client
c = cdsapi.Client()

# Define the list of variable, model and time range that you want to download
variable = "tas"
model_list = ["CNRM-CM6-1", "CanESM5", "EC-Earth3", "GFDL-ESM4", "IPSL-CM6A-LR", "MIROC-ES2L", "MPI-ESM1-2-LR",
              "MRI-ESM2-0"]
year_range = "1981-2010"


async def download_data(model):
    # Define the CDS API request
    request = {
        "variable": variable,
        "model": model,
        "experiment": "historical",
        "frequency": "day",
        "year": year_range,
        "format": "netcdf"
    }

    # Download the data using the CDS API client
    await c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': '2m_temperature',
        'year': '2017',
        'month': '01',
        'day': '01',
        'time': '12:00',
        'format': 'netcdf'
    },
    'download.zip')



async def main():
    tasks = []
    for model in model_list:
        tasks.append(asyncio.ensure_future(download_data(model)))
    await asyncio.gather(*tasks)


if __name__ == "__main__":
    loop = asyncio.get_event_loop()
    loop.run_until_complete(main())
