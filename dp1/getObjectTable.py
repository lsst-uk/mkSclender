import os
import re
import numpy as np
import pyvo
import requests
tap_url = 'https://data.lsst.cloud/api/tap'
token = 'gt-XXXXXXXXXXXXXXXXXXX' 
cred = pyvo.auth.CredentialStore()
cred.set_password("x-oauth-basic", token)
tap = pyvo.dal.TAPService(tap_url, session=cred.get("ivo://ivoa.net/sso#BasicAA"))

query = "SELECT count(*) " \
        "from dp1.Object"

results = tap.run_sync(query)
print(results.to_table())

idFile="dp1_objectIds.csv"
if not os.path.exists(idFile):
    query="SELECT objectId " \
          "from dp1.Object"
#    job = tap.run_async(query, output_format="csv")
    job = tap.submit_job(query, output_format="csv")
    job.run()
    if job.phase not in ("COMPLETED", "ERROR", "ABORTED"):
        job.wait()


    results = job.fetch_result()
    
    # Convert to Astropy Table and save CSV
    table = results.to_table()
    table.write(idFile, format="csv")
        
  

import pandas as pd


# Step 1: Read CSV (assume column name is 'objectId')
df = pd.read_csv(idFile)

# Step 2: Sort and deduplicate
object_ids = sorted(df["objectId"].astype(int).unique())

# Step 3: Split into 20 equal chunks by count
splits = np.array_split(object_ids, 50)
print(splits[0])
# Step 4: Build queries
queries = []
for split in splits:
    A, B = split[0], split[-1]   # range for this chunk
    query = f"select * from dp1.Object where objectId between {A} and {B};"
    queries.append(query)
print(queries[0])
print(queries[1])

loop=0
for q in queries:
    loop=loop+1
    print(q)
    resFile="dp1_objectIds_"+str(loop)+".csv"
    if not os.path.exists(resFile):
        job = tap.submit_job(q, output_format="csv")
        job.run()
        if job.phase not in ("COMPLETED", "ERROR", "ABORTED"):
            job.wait()


        results = job.fetch_result()
        table = results.to_table()
        table.write(resFile, format="csv")

