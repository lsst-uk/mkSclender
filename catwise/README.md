# catWISE - generating skinny tables for cross-match input

## Preparation

The input for the flagging was created running a bulk outgest on the catwise_2020 table in the WISE database held in the WFAU archives
The SQL is shown in [catwiseOutgest.sql](catwiseOutgest.sql)

## Skinny table flagging

The flagging described in [catwise_columns.pdf](catwise_columns.pdf) was then implemneted via [catwise3311.py](catwise311.py). The qserv partitioner was then used to chunk this
output.