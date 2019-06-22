from sklearn.decomposition import IncrementalPCA
import numpy as np
import csv
import pandas as pd
import sys

mode = int(sys.argv[1])

if mode == 0:
   path = "/dfs/scratch0/dataset/20161030-Boeing/data/WD673_rt1.txt"
   outfile = "boeing_pca.out"
elif mode== 1:
   path = "boeing_medium.out"
   outfile = "boeing_medium_result.out"
else:
   path = "boeing.out"
   outfile = "test.out"

BATCH_SIZE = 50000

sampleEvery = 10

def transformChunk(c):
    c.fillna(inplace=True, method="ffill")
    c.fillna(inplace=True, value=0)

dims = 13
pca = IncrementalPCA(n_components=dims)
skipRow = lambda x: x % sampleEvery != 0
useCol = lambda x: x != 't' and x != 'flight_index' 
reader = pd.read_csv(path, sep='\t', chunksize=BATCH_SIZE, skiprows=skipRow, usecols=useCol)
for i, chunk in enumerate(reader):
    transformChunk(chunk)
    print("fitting", i)
    pca.partial_fit(chunk)


mean = pca.mean_
stdv = np.sqrt(pca.var_)
stdv[stdv==0] = 1.0

transformed = None
reader = pd.read_csv(path, sep='\t', chunksize=BATCH_SIZE, skiprows=skipRow, usecols=useCol)
for i, chunk in enumerate(reader):
    print("transforming", i)
    transformChunk(chunk)
    chunk = pca.transform((chunk.values-mean)/stdv)
    if transformed is None: transformed = chunk
    else: transformed = np.vstack((transformed, chunk))


np.savetxt(outfile, transformed, delimiter=",")
print(transformed.shape)
