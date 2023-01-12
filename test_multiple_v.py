import cdsapi

c = cdsapi.Client()

c.retrieve(
    'projections-cmip6',
    {
        'format': 'zip',
        'temporal_resolution': 'monthly',
        'year': [
            '2015'
        ],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'experiment': 'ssp5_8_5',
        'variable': ['near_surface_air_temperature'],
        'model': 'access_cm2',
    },
    'download.zip')