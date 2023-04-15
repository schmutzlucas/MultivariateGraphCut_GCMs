from concurrent.futures import ThreadPoolExecutor
import cdsapi
import os
import calendar

dir = r'/mnt/w/LucasSchmutz/MultivariateGraphCut_GCMs/data/ERA5'
#dir = r'w:\LucasSchmutz\mHM\volta_data\raw_data\era5'
os.chdir(dir)

c = cdsapi.Client()

# Define a function to call the CDS API and retrieve data
def cds_api_call(variable, year, month, days, name):
    c.retrieve(
        'reanalysis-era5-single-levels',  # Dataset ID
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',  # Output format
            'year': year,  # Years to retrieve data for
            'month': month,  # Months to retrieve data for
            'day': days,  # Days to retrieve data for
            'time': time,
            'variable': variable,  # Variable name
        },
    name
    )

variables = ['total_precipitation']

years_int = list(range(1975,2024))
years = list(map(str, years_int))

months_int = list(range(1,13))
months = list(map(str, months_int))

month_31d = ['1','3','5','7','8','10','12']

time = ['00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00','08:00','09:00','10:00','11:00', 
        '12:00','13:00','14:00','15:00','16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00']

# Define the main function to run the data retrieval
def main():
    # create a thread pool with n worker threads
    with ThreadPoolExecutor(max_workers=10) as exe:
        z = 0
        for variable in variables:
            # Create a directory to store the downloaded files for each
            # model
            os.chdir(dir)
            var_dir = os.path.abspath(f"{variable}")
            if not os.path.exists(var_dir):
                os.makedirs(var_dir)
                print(f"Created directory: {var_dir}")
            os.chdir(var_dir)
            # Retrieve data for each variable
            for year in years:
                os.chdir(var_dir)
                for month in months:
                    if (month not in month_31d) and (month != '2'):
                        days_int = list(range(1,31))
                        days = list(map(str,days_int))
                    elif month in month_31d:
                        days_int = list(range(1,32))
                        days = list(map(str,days_int))
                    elif (month == '2') and (calendar.isleap(int(year)) == True):
                        days_int = list(range(1,30))
                        days = list(map(str,days_int))
                    elif (month == '2') and (calendar.isleap(int(year)) == False):
                        days_int = list(range(1,29))
                        days = list(map(str,days_int))
                    else:
                        print('ERROR: number of days in month not assignable')
                        exit
                    # Set the output filename for the downloaded file
                    name = os.path.abspath(f'{dir}\{variable}\era5_{variable}_{year}_{month}.nc')
                    # Check if the file already exists before downloading it
                    if not os.path.exists(name):
                        z += 1
                        try:
                            # Submit a job to the thread pool to download the
                            # data
                            var = [exe.submit(cds_api_call, variable, year,
                                                month, days, name)]
                            #print(f'{variable}_{year}_{month}_{max(map(int,days))}')
                        finally:
                            pass
                    else:
                        print(f'File era5_{variable}_{year}_{month}.nc already exists!')
    print('Finished!')

# Call the main function to start the data retrieval
main()

