import cdsapi
import asyncio

c = cdsapi.Client()


def cds_api_call(experiment,
                 model,
                 year, month,
                 variable,
                 name):
    c.retrieve(
        'projections-cmip6',
        {
            'format': 'zip',
            'temporal_resolution': 'monthly',
            'year': year,
            'month': month,
            'experiment': experiment,
            'variable': variable,
            'model': model,

        },
        name)


# Variables
month_int = list(range(1, 13))
month = list(map(str, month_int))
year_int = list(range(1950, 2015))
year = list(map(str, year_int))
print(year)

model_list = ['GFDL-ESM4',
              'GISS-E2-1-G',
              'FIO-ESM-2-0',
              'FGOALS-f3-L',
              'MPI-ESM1-2-HR',
              'CAMS-CSM1-0',
              'ACCESS-CM2',
              'AWI-CM-1-1-MR',
              'NM-CM5-0',
              'MPI-ESM1-2-LR',
              'MIROC6',
              'ACCESS-ESM1-5']

variable_list = ['tas', 'pr', 'tasmax', 'tasmin', 'huss']
# variable_list = ['pr', 'tas', 'huss']

experiment = ['historical']


def main():
    for n in experiment:
        for i in model_list:
            for j in variable_list:
                name = j + '_Amon_' + i + '_' + n + '.zip'
                try:
                    cds_api_call(n, i, year, month, j, name)
                finally:
                    pass


main()
