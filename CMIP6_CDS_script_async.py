import cdsapi
import asyncio

c = cdsapi.Client()


async def cds_api_call(experiment,
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
        name
    )


# Variables
month_int = list(range(1, 13))
month = list(map(str, month_int))
year_int = list(range(1950, 2015))
year = list(map(str, year_int))
print(year)

model_list = [  'AWI-CM-1-1-MR',
                'BCC-CSM2-MR',
                'CAMS-CSM1-0',
                'CanESM5',
                'CESM2-WACCM',
                'CESM2',
                'CIESM',
                'CMCC-CM2-SR5',
                'EC-Earth3-Veg',
                'EC-Earth3',
                'FGOALS-f3-L',
                'FGOALS-g3',
                'FIO-ESM-2-0',
                'GFDL-CM4',
                'GFDL-ESM4',
                'INM-CM4-8',
                'INM-CM5-0',
                'IPSL-CM6A-LR',
                'KACE-1-0-G',
                'KIOST-ESM',
                'MIROC6',
                'MPI-ESM1-2-HR',
                'MPI-ESM1-2-LR',
                'MRI-ESM2-0',
                'NESM3',
                'NorESM2-LM',
                'NorESM2-MM']

# for daily data
# variable_list = ['tas', 'pr', 'tasmax', 'tasmin', 'huss','uas', 'vas']

# for monthly data
variable_list = ['pr', 'tas', 'huss', 'uas', 'vas']

experiment = ['historical']


async def main():
    z = 0
    tasks = []
    for n in experiment:
        for i in model_list:
            for j in variable_list:
                name = 'download/' + j + '_Amon_' + i + '_' + n + '.zip'
                z += 1
                try:
                    tasks.append(cds_api_call(n, i, year, month, j, name))
                finally:
                    pass
    await asyncio.gather(*tasks)


asyncio.run(main())
