# casc
Solver for Context-Aware Segmentation and Clustering for Motif Discovery in Noisy Time Series Data

## Instructions

From the main directory import the file CASC_solver `from CASC_solver import CASCSolver`

Then create a solver object. The solver has the following options:
```
solver = CASCSolver(
  window_size,
  number_of_clusters,
  lambda_parameter,
  beta,
  threshold, # convergence threshold
  gamma,
  input_file, # input data file
  num_proc, # number of processes running in parallel as workers
  maxMotifs, # cap number of motifs
  motifReq, # minimum number of motifs
  maxIters, # number of iterations to run (None if until convergence)
  )
```
Then use the solver to run CASC:
```
(cluster_assignment, cluster_MRFs, motifs, motifRanked, bic, runtime) = solver.PerformFullCASC(
        initialClusteredPoints, # the initial clustered points if you want to start with a pre-assignment
        useMotif # whether to use motifs (if false then just performs TICC until convergence)
```

The input data file should be a csv with one line per time step and each line having the sensor values for that step. This file can be PCA'd down if necessary. The output files will be a `cluster` -> a list of primary cluster labels given per time step, `cluster_MRFs` -> the inverse covariance matrices learned, `motifs` -> the motifs found as well as their identified instances, `motifsRanked` -> the scores for each motif.

## Directory Structure
To run paper code, extract scripts in the scripts folder into the main directory. Synthetic data was generated via the generate_synthetic.py script. The synthetic experiment scripts are in the scripts/synthetic folder. Case study scripts are in scripts/case_studies. The plotting and graphing scripts are in scripts/aggregation_and_plotting. These scripts assume that the zipped data in synthetic has been unzipped and placed into the main folder.
