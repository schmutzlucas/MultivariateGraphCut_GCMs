import cdsapi
import asyncio


def background(f):
    def wrapped(*args):
        return asyncio.get_event_loop().run_in_executor(None, f, *args)
    return wrapped


@background
def cds_api_call(experiment,
                 model,
                 year, month,
                 variable,
                 name):
    c.retrieve(
        'projections-cmip6',
        {
            'temporal_resolution': 'monthly',
            'experiment': experiment,
            'variable': variable,
            'model': model,
            'format': 'zip',
            'year': year,
            'month': month,
        },
        name)


model_list = ['access_esm1_5']

variable_list = [
    'tas'
]
# variable_list = ['pr', 'tas', 'huss']

scenario_list = [
    'historical',
    'ssp5_8_5'
]

for n in scenario_list:
    for i in model_list:
        for j in variable_list:
            name = j + '_Amon_' + i + '_' + n + '.zip'
            try:
                cds_api_call(n, i, j, name)
            finally:
                pass