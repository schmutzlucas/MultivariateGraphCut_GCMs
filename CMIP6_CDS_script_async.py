from concurrent.futures import ThreadPoolExecutor

import cdsapi

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
        name
    )


with open('model_list.txt', 'r') as file:
    model_list = file.read().splitlines()


# for daily data
# variable_list = ['tas', 'pr', 'tasmax', 'tasmin', 'huss','uas', 'vas']

# for monthly data
variable_list = ['pr', 'tas', 'huss', 'uas', 'vas', 'tasmax', 'tasmin', 'psl', 'prw']

experiment = ['historical', 'ssp5_8_5']


def main():
    # create a thread pool with n worker threads
    with ThreadPoolExecutor(max_workers=10) as exe:
        z = 0
        for n in experiment:
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
                for j in variable_list:
                    name = 'download/' + j + '_Amon_' + i + '_' + n + '.zip'
                    z += 1
                    try:
                        var = [exe.submit(cds_api_call, n, i, year, month, j, name)]
                    finally:
                        pass


main()
