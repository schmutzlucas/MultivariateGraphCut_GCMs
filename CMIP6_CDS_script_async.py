import cdsapi
import asyncio

c = cdsapi.Client()


async def cds_api_call(experiment,
                       model,
                       year, month,
                       variable,
                       name):
    await c.retrieve(
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

model_list = ['ACCESS-CM2',
              'ACCESS-ESM1-5',
              'AWI-CM-1-1-MR',
              'AWI-ESM-1-1-LR',
              'BCC-CSM2-MR',
              'BCC-ESM1',
              'CAMS-CSM1-0',
              'CanESM5',
              'CanESM5-CanOE',
              'CESM2 ',
              'CESM2-FV2 ',
              'CESM2-WACCM ',
              'CESM2-WACCM-FV2 ',
              'CIESM',
              'CMCC-CM2-HR4',
              'CMCC-CM2-SR5',
              'CMCC-ESM2',
              'CNRM-CM6-1',
              'CNRM-CM6-1-HR',
              'CNRM-ESM2-1',
              'E3SM-1-0 ',
              'E3SM-1-1 ',
              'E3SM-1-1-ECA ',
              'EC-Earth3',
              'EC-Earth3-AerChem',
              'EC-Earth3-CC',
              'EC-Earth3-Veg',
              'EC-Earth3-Veg-LR',
              'FGOALS-f3-L',
              'FGOALS-g3',
              'FIO-ESM-2-0',
              'GFDL-ESM4 ',
              'GISS-E2-1-G ',
              'GISS-E2-1-H ',
              'HadGEM3-GC31-LL',
              'HadGEM3-GC31-MM',
              'IITM-ESM',
              'INM-CM4-8',
              'INM-CM5-0',
              'IPSL-CM5A2-INCA',
              'IPSL-CM6A-LR',
              'KACE-1-0-G',
              'KIOST-ESM',
              'MCM-UA-1-0 ',
              'MIROC6',
              'MIROC-ES2H',
              'MIROC-ES2L',
              'MPI-ESM-1-2-HAM',
              'MPI-ESM1-2-HR',
              'MPI-ESM1-2-LR',
              'MRI-ESM2-0',
              'NESM3',
              'NorCPM1',
              'NorESM2-LM',
              'NorESM2-MM',
              'SAM0-UNICON',
              'TaiESM1',
              'UKESM1-0-LL']
model_list = ['IITM-ESM']

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
