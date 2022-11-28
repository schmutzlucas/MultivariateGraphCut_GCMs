import cdsapi
import asyncio


def background(f):
    def wrapped(*args):
        return asyncio.get_event_loop().run_in_executor(None, f, *args)

    return wrapped


@background
def cds_api_call(n, i, j, name):
    cdsapi.Client().retrieve(
        'projections-cmip6',
        {
            'temporal_resolution': 'Monthly',
            'experiment': n,
            'level': 'single_levels',
            'variable': j,
            'model': i,
            'format': 'zip'
        },
        name)

@background
def cds_api_call(n, i, j, name):
    c.retrieve(
        'projections-cmip6',
        {
            'temporal_resolution': 'monthly',
            'experiment': 'ssp5_8_5',
            'variable': 'near_surface_air_temperature',
            'model': 'access_cm2',
            'format': 'zip',
            'year': [
                '2015', '2016', '2017',
                '2018', '2019', '2020',
                '2021', '2022', '2023',
                '2024', '2025', '2026',
                '2027', '2028', '2029',
                '2030', '2031', '2032',
                '2033', '2034', '2035',
                '2036', '2037', '2038',
                '2039', '2040', '2041',
                '2042', '2043', '2044',
                '2045', '2046', '2047',
                '2048', '2049', '2050',
                '2051', '2052', '2053',
                '2054', '2055', '2056',
                '2057', '2058', '2059',
                '2060', '2061', '2062',
                '2063', '2064', '2065',
                '2066', '2067', '2068',
                '2069', '2070', '2071',
                '2072', '2073', '2074',
                '2075', '2076', '2077',
                '2078', '2079', '2080',
                '2081', '2082', '2083',
                '2084', '2085', '2086',
                '2087', '2088', '2089',
                '2090', '2091', '2092',
                '2093', '2094', '2095',
                '2096', '2097', '2098',
                '2099', '2100',
            ],
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
        },
        'download.zip')


# model_list = [
#     'gfdl_esm4',
#     'giss_e2_1_h',
#     'fio_esm_2_0',
#     'fgoals_f3_l',
#     'mpi_esm1_2_hr',
#     'cams_csm1_0',
#     'access_cm2',
#     'awi_cm_1_1_mr',
#     'inm_cm5_0',
#     'mpi_esm1_2_lr',
#     'miroc6',
#     'access_esm1_5',
#     'bcc_csm2_mr',
#     'kace_1_0_g',
#     'mcm_ua_1_0'
# ]
model_list = ['access_esm1_5']

# variable_list = [
#     'near_surface_air_temperature',
#     'daily_maximum_near_surface_air_temperature',
#     'daily_minimum_near_surface_air_temperature',
#     'precipitation',
#     'sea_level_pressure',
#     'near_surface_relative_humidity',
#     'near_surface_specific_humidity',
#     'near_surface_wind_speed'
# ]
# variable_list = [
#     'tas',
#     'near_surface_wind_speed'
# ]

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

print('download finished')
