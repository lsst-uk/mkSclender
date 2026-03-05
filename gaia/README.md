# GAIADR3 - generating skinny tables for cross-match input

## Preperation

The input for the flagging was created running a bulk outgest on the gaia_source table in the GAIADR3 database held in the WFAU archives
The SQL is shown in [gaiaDR3Outgest.sql](gaiaDR3Outgest.sql)

## Skinny table flagging

The flagging described in [gaia_dr3.pdf](gaia_dr3.pdf) was then implemneted via [gaiaDR3311.py](gaiaDR3311.py). The qserv partitioner was then used to chunk this
output.