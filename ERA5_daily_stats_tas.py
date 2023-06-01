from concurrent.futures import ThreadPoolExecutor
import cdsapi
import os
import calendar

dir = dir = r'/mnt/w/LucasSchmutz/MultivariateGraphCut_GCMs/data/ERA5/raw'
os.chdir(dir)

c = cdsapi.Client()

def cds_api_call_new_api(variable, year, month, name):
    result = c.service(
        "tool.toolbox.orchestrator.run_workflow.user-apps.app-c3s-daily-era5-statistics",
        kwargs={
            "dataset": "reanalysis-era5-single-levels",
            "product_type": "reanalysis",
            "variable_e5sl": variable,
            "statistic": "daily_mean",
            "year_e5sl": year,
            "month": month,
            "time_zone": "UTC+00:00",
            "frequency": "1-hourly",
            "grid_e5": "3.0/3.0",
            "pressure_level_e5sl": "-",
            "area": {"lat": [-90, 90], "lon": [-180, 180]},
        },
        workflow_name="application",
    )
    c.download(result)


variables = ['2m_temperature']

years_int = list(range(1975, 2015))
years = list(map(str, years_int))

months_int = list(range(1, 13))
months = list(map(str, months_int))


def main():
    with ThreadPoolExecutor(max_workers=10) as exe:
        for variable in variables:
            os.chdir(dir)
            var_dir = os.path.abspath(f"{variable}")
            if not os.path.exists(var_dir):
                os.makedirs(var_dir)
                print(f"Created directory: {var_dir}")
            os.chdir(var_dir)

            for year in years:
                for month in months:
                    if variable == '2m_temperature':
                        var_short = 'tas'
                        model = 'ERA5'
                    elif variable == 'total_precipitation':
                        var_short = 'pr'
                        model = 'ERA5'
                    else:
                        print(f"ERROR: Invalid variable {variable}")
                        exit()

                    name = os.path.abspath(f'{dir}/{variable}/{var_short}_{model}_daily_mean_{year}{month}_regridded.nc')


                    if not os.path.exists(name):
                        exe.submit(cds_api_call_new_api, variable, year, month, name)
                    else:
                        print(f'File era5_{variable}_{year}_{month}_daily_mean.nc already exists!')

    print('Finished!')

main()
