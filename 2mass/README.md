# 2MASS - generating skinny tables for cross-match input

## Preparation

The input for the flagging was created running a bulk outgest on the twomass_psc and twomass_xsc tables in the TWOMASS database held in the WFAU archives
The SQL is shown in [2massOutgest.sql](2massOutgest.sql). Only rows with no entry in twomass_xsc were outgested from twomass_psc.

## Skinny table flagging

The flagging described in [2mass.pdf](2mass.pdf) was then implemneted via [psc311.py](psc311.py) and [xsc311.py](xsc311.py). The qserv partitioner was then used to chunk this
output.