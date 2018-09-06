from CASC_solver import CASCSolver
import numpy as np
import sys
import pickle
import os

CLUSTER_NUMBER = 10

def pickleObject(fname, data):
    f = open(fname, "wb")
    pickle.dump(data, f)
    f.close()

def runCASC(inputName, outputDir):
    solver = CASCSolver(window_size=1, number_of_clusters=10,
        lambda_parameter=1e-3, beta=25, threshold=2e-5,
        gamma=0.8, input_file=inputName, num_proc=10, maxMotifs=50,
        motifReq=10, maxIters=5)
    (cluster_assignment, cluster_MRFs, motifs, motifRanked, bic, _) = solver.PerformFullCASC(useMotif=True)
    solver.CleanUp()
    motifFile = "%smotifs.pkl" % outputDir
    pickleObject(motifFile, motifs)
    motifRankedFile = "%smotifRanked.pkl" % outputDir
    pickleObject(motifRankedFile, motifRanked)
    outputName = None
    if outputDir is not None:
        outputName = "%sassign.out" % outputDir
        np.savetxt(outputName, cluster_assignment, fmt='%d')
    return outputName, bic


if __name__ == "__main__":
    # mode of 1 to skip old assign
    assert len(sys.argv) == 3
    input_fname = sys.argv[1]
    output_fdir = sys.argv[2]
    runCASC(input_fname, output_fdir)
