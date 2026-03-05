# GAIADR3 - generating skinny tables for cross-match input

## Preperation

The input for the flagging was created running a bulk outgest on the gaia_source table in the GAIADR3 database held in the WFAU archives
The SQL is shown in [gaiaDR3Outgest.sql](gaiaDR3Outgest.sql)

As DP1 is small the full table was outgested into 50 CSV files (select * from Objecct where objectId between x and Y).
[extractColsToCSV.py](extractColsToCSV.py) was then used to pick out the required columns and combine into a single CSV file.

## Skinny table flagging

The flagging described in [lsst_dp1.pdf](lsst_dp1.pdf) was then implemneted via [lsst311.py](lsst311.py)