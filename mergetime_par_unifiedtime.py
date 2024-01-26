import cftime
import pandas as pd
import os
import xarray as xr
import numpy as np
from datetime import datetime
from multiprocessing import Pool
from cftime import DatetimeProlepticGregorian

def get_time_constraint(start, end):
    start_date = pd.Timestamp(start).to_pydatetime()
    end_date = pd.Timestamp(end).to_pydatetime()

    time_constraint = iris.Constraint(time=lambda cell: start_date <= cell.point <= end_date)
    return time_constraint

root_dir = 'Y:\LucasSchmutz\MultivariateGraphCut_GCMs\download_day_unzip'
merged_dir = 'Y:\LucasSchmutz\MultivariateGraphCut_GCMs\download_day_merged_multithread'

def preprocess(ds):
    ds = ds.copy()
    time_units = "days since 1850-01-01"
    calendar = "proleptic_gregorian"

    decoded_time = xr.coding.times.decode_cf_datetime(ds['time'],
                                                      units=time_units,
                                                      calendar=calendar)
    ds['time'] = xr.DataArray(
        [cftime.datetime(t.year, t.month, t.day) for t in decoded_time],
        dims='time')
    return ds

def merge_model_dir(model_dir):
    for var_exp_dir in os.listdir(os.path.join(root_dir, model_dir)):
        if not os.path.isdir(os.path.join(root_dir, model_dir, var_exp_dir)):
            continue

        output_dir = os.path.join(merged_dir, model_dir, var_exp_dir)
        os.makedirs(output_dir, exist_ok=True)

        input_files = []
        for filename in os.listdir(
                os.path.join(root_dir, model_dir, var_exp_dir)):
            if filename.endswith(".nc"):
                input_file_path = os.path.join(root_dir, model_dir, var_exp_dir,
                                               filename)
                input_files.append(input_file_path)

        print(f"Merging {len(input_files)} files:")
        print('\n'.join(input_files))
        now = datetime.now()
        print("Current time:", now)

        try:
            ds = xr.open_mfdataset(input_files, combine='nested',
                                   concat_dim='time', preprocess=preprocess)

            # Apply time constraint here
            start_date = DatetimeProlepticGregorian(
                pd.Timestamp(start).to_pydatetime())
            end_date = DatetimeProlepticGregorian(
                pd.Timestamp(end).to_pydatetime())
            time_constraint = iris.Constraint(
                time=lambda cell: start_date <= cell.point <= end_date)

            ds = ds.sel(time=time_constraint)

            start_date = ds.time.values[0].strftime('%Y%m%d')
            end_date = ds.time.values[-1].strftime('%Y%m%d')

            input_filename = os.path.basename(input_files[0])
            var = input_filename.split("_")[0]
            model = input_filename.split("_")[2]
            exp = input_filename.split("_")[3]
            ens = input_filename.split("_")[4]
            version = input_filename.split("_")[-1].split(".")[0]
            output_filename = f"{var}_{model}_{exp}_{ens}_{start_date}-{end_date}_{version}.nc"

            output_path = os.path.join(output_dir, output_filename)
            print(output_path)
            if os.path.exists(output_path):
                print(
                    f"Output file {output_filename} already exists, skipping...")
            else:
                ds.to_netcdf(output_path)

        except Exception as e:
            current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            with open("failed_merge.log", "a") as f:
                f.write(
                    f"{current_time} - An error occurred while processing {input_files}: {e}\n")
            print(f"An error occurred while processing {input_files}: {e}")


if __name__ == '__main__':
    model_dirs = [d for d in os.listdir(root_dir) if
                  os.path.isdir(os.path.join(root_dir, d))]
    num_processes = 8
    pool = Pool(num_processes)
    pool.map(merge_model_dir, model_dirs)
    pool.close()
    pool.join
    print("Script complete!")

