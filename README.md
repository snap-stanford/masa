# MASA
Solver for MASA: Motif-Aware State Assignment (previously called CASC)

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
The code from the paper is in the directory `paper_code`. To run a script, put the script in the main directory.

### Synthetic Data Code
The scripts for the synthetic experiments are in `paper_code/scripts/synthetic`. `baseline.py` contains the script for running the baselines, while `synthetic.py` contains the script for running MASA.

The synthetic data can be found in `ordered_synthetic.zip`. You need to unzip that file and put it in the main directory. The script that was used to generate that data is found in `paper_code/generateDatasets/generate_synthetic.py.`

### Cycling data
The cycling data can be found in cycling.zip. The script to create the cycling data is in `scripts/cycling/create_cycling_dataset.py` and the script to run the cycling data with MASA is in `scripts/cycling.py`. The actual cycling data is in `cycling.zip`.

### Case Studies
Unfortunately we cannot release the datasets for the automobile and airplane data. The scripts that were used to run MASA on this data can be found in `paper_code/scripts/runCaseStudy.py`.

### Aggregation and plotting
Aggregation and plotting scripts can be found in `paper_code/scripts/aggregation_and_plotting`