from concurrent.futures import ThreadPoolExecutor
import cdsapi
import os

# Create a CDS API client object
c = cdsapi.Client()


# Define a function to call the CDS API and retrieve data
def cds_api_call(experiment, model, year, month, variable, name):
    c.retrieve(
        'projections-cmip6',  # Dataset ID
        {
            'format': 'zip',  # Output format
            'temporal_resolution': 'monthly',  # Temporal resolution of the data
            'year': year,  # Years to retrieve data for
            'month': month,  # Months to retrieve data for
            'experiment': experiment,  # Experiment name
            'variable': variable,  # Variable name
            'model': model,  # Model name
        },
        name  # Output filename
    )


# Read the list of models from a file
with open('../model_list.txt', 'r') as file:
    model_list = file.read().splitlines()

# Define the list of variables to retrieve
variable_list = ['pr', 'tas', 'huss', 'uas', 'vas',
                 'tasmax', 'tasmin', 'psl', 'prw']

# Define the list of experiments to retrieve data for
experiment = ['historical', 'ssp5_8_5']


# Define the main function to run the data retrieval
def main():
    # create a thread pool with n worker threads
    with ThreadPoolExecutor(max_workers=10) as exe:
        z = 0
        for n in experiment:
            # Set the year and month ranges depending on the experiment
            if n == 'historical':
                month_int = list(range(1, 13))
                month = list(map(str, month_int))
                year_int = list(range(1950, 2015))
                year = list(map(str, year_int))
            elif n == 'ssp5_8_5':
                month_int = list(range(1, 13))
                month = list(map(str, month_int))
                year_int = list(range(2015, 2101))
                year = list(map(str, year_int))
            for i in model_list:
                # Create a directory to store the downloaded files for each
                # model
                directory = os.path.abspath(f"download/{i}")
                if not os.path.exists(directory):
                    os.makedirs(directory)
                    print(f"Created directory: {directory}")
                # Retrieve data for each variable
                for j in variable_list:
                    # Set the output filename for the downloaded file
                    name = f"download/{i}/{j}_Amon_{i}_{n}.zip"
                    # Check if the file already exists before downloading it
                    if not os.path.exists(name):
                        z += 1
                        try:
                            # Submit a job to the thread pool to download the
                            # data
                            var = [exe.submit(cds_api_call, n, i, year,
                                              month, j, name)]
                        finally:
                            pass
    print('finished')


# Call the main function to start the data retrieval
main()
