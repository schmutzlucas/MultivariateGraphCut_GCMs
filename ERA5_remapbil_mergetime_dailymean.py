import xarray as xr
import numpy as np
from scipy.interpolate import griddata
import glob

# Define the grid from my_grid.txt
# Replace the following two lines with the content of your my_grid.txt file
lons = np.linspace(min_lon, max_lon, num_lon_points)
lats = np.linspace(min_lat, max_lat, num_lat_points)

# Read and merge the data
files = glob.glob('ERA5*temp*.nc')
ds = xr.open_mfdataset(files, combine='by_coords')

# Remap the data to the new grid
lons_new, lats_new = np.meshgrid(lons, lats)
data_remapped = ds['t2m'].groupby('time').apply(
    lambda x: griddata(
        (ds.longitude.values.ravel(), ds.latitude.values.ravel()),
        x.values.ravel(),
        (lons_new, lats_new),
        method='linear'
    )
)

# Create a new xarray Dataset with the remapped data
ds_remapped = xr.Dataset(
    {
        't2m': (('time', 'latitude', 'longitude'), data_remapped)
    },
    coords={
        'time': ds.time,
        'latitude': (('latitude',), lats),
        'longitude': (('longitude',), lons)
    }
)

# Compute daily means
ds_daily = ds_remapped.resample(time='D').mean()

# Save the output to a NetCDF file
ds_daily.to_netcdf('tas_ERA5_1975-2024_daily_remapped.nc')
