import os
import xarray as xr
import numpy as np
from datetime import datetime
from multiprocessing import Pool
from cftime import DatetimeProlepticGregorian

root_dir = 'Y:\LucasSchmutz\MultivariateGraphCut_GCMs\download_day_unzip'
merged_dir = 'Y:\LucasSchmutz\MultivariateGraphCut_GCMs\download_day_merged_multithread'


def preprocess(ds):
    ds = ds.copy()
    time_units = "days since 1850-01-01"
    calendar = "proleptic_gregorian"

    ds['time'] = xr.coding.times.decode_cf_datetime(ds['time'],
                                                    units=time_units,
                                                    calendar=calendar)
    ds['time'] = ds.indexes['time'].normalize()

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
            start_date = np.datetime_as_string(ds.time.values[0],
                                               unit='D').replace('-', '')
            end_date = np.datetime_as_string(ds.time.values[-1],
                                             unit='D').replace('-', '')

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
            with open("failed_merge.log", "a") as f:
                f.write(
                    f"An error occurred while processing {input_files}: {e}\n")
            print(f"An error occurred while processing {input_files}: {e}")


if __name__ == '__main__':
    model_dirs = [d for d in os.listdir(root_dir) if
                  os.path.isdir(os.path.join(root_dir, d))]
    num_processes = 8
    pool = Pool(num_processes)
    pool.map(merge_model_dir, model_dirs)
    pool.close()
    pool.join()
    print("Script complete!")
