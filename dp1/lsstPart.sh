sph-partition --config-file=/disk85/hips/DP1/lsstPart.json  --in.path=dp1Good.csv --mr.num-workers=2  --mr.pool-size=4096 --mr.block-size=16  --verbose --out.dir=/disk85/hips/DP1/chunks/ >> part.out

sph-partition --config-file=/disk85/hips/DP1/lsstPart.json  --in.path=dp1Bad.csv --mr.num-workers=2  --mr.pool-size=4096 --mr.block-size=16  --verbose --out.dir=/disk85/hips/DP1/rejectChunks/ >> part.out



#sph-partition --config-file=/disk85/hips/2MASS/pscPart.json  --in.path=pscGood.csv --mr.num-workers=2  --mr.pool-size=4096 --mr.block-size=16  --verbose --out.dir=/disk85/hips/2MASS/chunks/ >> part.out
