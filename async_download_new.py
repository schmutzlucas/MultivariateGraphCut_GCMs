import cdsapi
import asyncio

c = cdsapi.Client()


#
# def background(f):
#     def wrapped(*args):
#         return asyncio.get_event_loop().run_in_executor(None, f, *args)
#
#     return wrapped
#
#
# @background
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
        name)


# Variables
month = list(range(1, 13))
year = list(range(1950, 2015))

model_list = ['ACCESS-CM2']

variable_list = ['tas', 'pr']
# variable_list = ['pr', 'tas', 'huss']

experiment = ['historical', 'ssp5_8_5']

month = list(range(1, 13))
year = list(range(1950, 2015))

for n in experiment:
    for i in model_list:
        for j in variable_list:
            name = j + '_Amon_' + i + '_' + n + '.zip'
            try:
                cds_api_call(n, i, year, month, j, name)
            finally:
                pass
