# DP1 - generating skinny tables for cross-match input

## Preperation

As DP1 was not available for ingestion into the UK-DAC, the DP1 Object table was first outgested using
[getObjectTable.py](getObjectTable.py) which made a series of queries
to the TAP service on the US RSP.

As DP1 is small the full table was outgested into 50 CSV files (select * from Objecct where objectId between x and Y).
[extractColsToCSV.py](extractColsToCSV.py) was then used to pick out the required columns and combine into a single CSV file.

## Skinny table flagging

The flagging described in [lsst_dp1.pdf](lsst_dp1.pdf) was then implemneted via [lsst311.py]