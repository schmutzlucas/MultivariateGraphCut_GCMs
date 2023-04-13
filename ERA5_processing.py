import xarray as xr
from multiprocessing import Pool
import glob
import os

vars = ['pre','temp']

def process_variable(var):
    print(f'Processing variable: {var}...')
    if var == 'pre':
        inDir  = 'X:\\LoicGerber\\data\\raw_data\\era5\\total_precipitation\\'
        outDir = 'X:\\LoicGerber\\data\\processed_data\\era5\\total_precipitation\\'
        if not os.path.exists(outDir):
            os.makedirs(os.path.join(outDir,'daily'))
            print(f'Output directory for {var} created')
    elif var == 'temp':
        # Define the directory where your netcdf files are located
        inDir  = 'X:\\LoicGerber\\data\\raw_data\\era5\\2m_temperature\\'
        outDir = 'X:\\LoicGerber\\data\\processed_data\\era5\\2m_temperature\\'
        if not os.path.exists(outDir):
            os.makedirs(os.path.join(outDir,'daily'))
            print(f'Output directory for {var} created')
    # Use glob to get a list of all the netcdf files in the directory
    files = sorted(glob.glob(inDir + '*.nc'))
    # Open the first file to get the variable name and coordinates
    ds0 = xr.open_dataset(files[0])
    var_name = list(ds0.data_vars.keys())[0]
    # Loop through each file and calculate the daily sum
    for file in files:
        ds = xr.open_dataset(file)
        name = file.replace('.nc','')
        name = name.replace(inDir,'')
        if var == 'pre':
            name = name + '_daily'
            daily = ds.resample(time='D').sum()*1000.0 # to have mm
            # Rename the variable to the original variable name
            daily = daily.rename({var_name: 'pre'})
            daily['pre'].attrs['units'] = 'mm'
            print(f'Saving {name}...')
            daily.to_netcdf(os.path.join(outDir,'daily',f'{name}.nc'))
        elif var == 'temp':
            name_tavg = name + '_daily_tavg'
            name_tmax = name + '_daily_tmax'
            name_tmin = name + '_daily_tmin'
            daily_tavg = ds.resample(time='D').mean()-273.15
            daily_tmax = ds.resample(time='D').max()-273.15
            daily_tmin = ds.resample(time='D').min()-273.15
            daily_tavg = daily_tavg.rename({var_name: 'tavg'})
            daily_tmax = daily_tmax.rename({var_name: 'tmax'})
            daily_tmin = daily_tmin.rename({var_name: 'tmin'})
            daily_tavg['tavg'].attrs['units'] = '°C'
            daily_tmax['tmax'].attrs['units'] = '°C'
            daily_tmin['tmin'].attrs['units'] = '°C'
            for data, file_name in zip([daily_tavg, daily_tmax, daily_tmin], [name_tavg, name_tmax, name_tmin]):
                print(f'Saving {file_name}...')
                data.to_netcdf(os.path.join(outDir, 'daily', f'{file_name}.nc'))
    if var == 'pre':
        # Merge and save the merged file
        print(f'Saving global {var} files...')
        with xr.open_mfdataset(glob.glob(os.path.join(outDir, 'daily', '*.nc')), combine='nested',concat_dim='time') as daily_sum:
            daily_sum.to_netcdf(os.path.join(outDir,f'{var}.nc'))
    elif var == 'temp':
        # Merge and save the merged file
        print(f'Saving global {var} files...')
        with xr.open_mfdataset(glob.glob(os.path.join(outDir, 'daily', '*tavg.nc')), combine='nested',concat_dim='time') as daily_tavg:
            daily_tavg.to_netcdf(os.path.join(outDir,f'{var}.nc'))
        with xr.open_mfdataset(glob.glob(os.path.join(outDir, 'daily', '*tmax.nc')), combine='nested',concat_dim='time') as daily_tmax:
            daily_tmax.to_netcdf(os.path.join(outDir,f'{var}.nc'))
        with xr.open_mfdataset(glob.glob(os.path.join(outDir, 'daily', '*tmin.nc')), combine='nested',concat_dim='time') as daily_tmin:
            daily_tmin.to_netcdf(os.path.join(outDir,f'{var}.nc'))
    print(f'Finished processing variable {var}')

if __name__ == '__main__':
    pool = Pool(processes=len(vars))
    pool.map(process_variable, vars)
    pool.close()
    pool.join()

