import os
from datetime import datetime, timedelta

# Function to get the filenames for a given date range
def generate_filenames(start_date, end_date):
    current_date = start_date
    while current_date <= end_date:
        yield f"tas_ERA5_{current_date.strftime('%Y%m01')}-{current_date.strftime('%Y%m31')}.nc"
        current_date += timedelta(days=31)
        current_date = current_date.replace(day=1)

# Function to check the existence of files in a directory
def check_files(directory, start_date, end_date):
    missing_files = []
    for filename in generate_filenames(start_date, end_date):
        if not os.path.isfile(os.path.join(directory, filename)):
            missing_files.append(filename)
    return missing_files

# Testing the function
directory = "data/ERA5/tas/"  # REPLACE WITH YOUR PATH
start_date = datetime.strptime("1960-01-01", "%Y-%m-%d")
end_date = datetime.strptime("2022-12-31", "%Y-%m-%d")

missing_files = check_files(directory, start_date, end_date)

if len(missing_files) == 0:
    print("All files are here.")
else:
    print("The following files are missing:")
    for file in missing_files:
        print(file)
