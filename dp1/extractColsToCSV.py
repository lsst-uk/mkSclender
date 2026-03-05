import os
import pandas as pd
from pprint import pprint

input_dir = "/disk85/hips/DP1"          # directory containing dp1_objectIds_*.csv
output_file = "dp1_combined_mags.csv"         # final output filename

base_cols = ["objectId","coord_ra", "coord_dec","coord_raErr","coord_decErr","coord_ra_dec_Cov"]

# The filters you want to expand over
filters = ["u", "g", "r", "i", "z", "y"]

# Start list with the base columns in-order
columns_to_use = base_cols.copy()

# Loop over filters and append mag + magerr columns for each
for f in filters:
    columns_to_use.append(f"{f}_cModelMag")
    columns_to_use.append(f"{f}_cModelMagErr")
    columns_to_use.append(f"{f}_cModelFlux")
    columns_to_use.append(f"{f}_cModelFluxErr")
    columns_to_use.append(f"{f}_cModel_flag")
    columns_to_use.append(f"{f}_psfMag")
    columns_to_use.append(f"{f}_psfMagErr")
    columns_to_use.append(f"{f}_psfFlux")
    columns_to_use.append(f"{f}_psfFluxErr")
    columns_to_use.append(f"{f}_psfFlux_flag")
    columns_to_use.append(f"{f}_pixelFlags_edge")
    columns_to_use.append(f"{f}_pixelFlags_offimage")
    columns_to_use.append(f"{f}_pixelFlags_saturated")
    columns_to_use.append(f"{f}_pixelFlags_bad")
    columns_to_use.append(f"{f}_pixelFlags_clipped")
    columns_to_use.append(f"{f}_pixelFlags_cr")
    columns_to_use.append(f"{f}_pixelFlags_sensor_edge")
    columns_to_use.append(f"{f}_pixelFlags_suspect")
    columns_to_use.append(f"{f}_epoch")    

pprint(columns_to_use)




# remove output file if it already exists (avoid accidental appends)
if os.path.exists(output_file):
    os.remove(output_file)

first_file = True

for fname in sorted(os.listdir(input_dir)):
    if fname.startswith("dp1_objectIds_") and fname.endswith(".csv"):
        fullpath = os.path.join(input_dir, fname)
        print(f"Processing: {fullpath}")

        # stream read in case files are large
        for chunk in pd.read_csv(fullpath, usecols=columns_to_use, chunksize=500000):
            chunk = chunk.reindex(columns=columns_to_use)
            chunk.to_csv(output_file,
                         index=False,
                         mode="a",
                         header=first_file)   # write header only once
            first_file = False
