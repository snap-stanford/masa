from TICC_solver import TICCSolver
import numpy as np
import sys
import pickle
import os

CLUSTER_NUMBER = 8

def dataset(mode, input_name, output_dir):
    beta = 100 
    number_of_clusters = CLUSTER_NUMBER
    if mode == 1:
        outputName = "%s/old/assign.out" % output_dir
    else:
        outputName, bic = runNonMotifTICC(input_name, output_dir, number_of_clusters, beta, None)
        print(bic)
    if mode == 1: runHyperParameterTests(input_name, output_dir, number_of_clusters, beta, outputName)


def runHyperParameterTests(inputName, outputDir, clusters, beta, oldAssignmentsName):
    gammas = [0.6, 0.8]
    motifReqs = 10
    for g in gammas:
        gammaDir = "%s/%s/" % (outputDir, g)
        makeDir(gammaDir)
        runTest(1, inputName, gammaDir, clusters,
                beta, g, motifReqs, oldAssignmentsName, 20)

def runBICTests(inputName, number_of_clusters):
    beta = [15, 25, 50, 75, 100, 500]
    bicBeta = []
    for b in beta:
        _, bic = runNonMotifTICC(inputName, None, number_of_clusters, b, None)
        bicBeta.append((bic, b))
    bicBeta.sort(reverse=True)
    print(bicBeta)


def runNonMotifTICC(inputName, outputDir, clusters, beta, oldAssignmentsName):
    if outputDir is not None: 
        oldDir = "%s/old/" % outputDir
        makeDir(oldDir)
        outputDir = oldDir
    return runTest(0, inputName, outputDir, clusters, beta, 1, 1, oldAssignmentsName, 40)

def pickleObject(fname, data):
    f = open(fname, "wb")
    pickle.dump(data, f)
    f.close()

def makeDir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def runTest(mode, inputName, outputDir, clusters, beta, gamma, motifReq, oldAssignmentsName, maxIters):
    print("TESTING %s" % (gamma))
    #maxIters used to be 30
    solver = TICCSolver(window_size=1, number_of_clusters=clusters, lambda_parameter=1e-3, beta=beta, threshold=2e-5,
                        gamma=gamma, input_file=inputName, num_proc=10, maxMotifs=50, motifReq=motifReq, maxIters=maxIters)
    old_assign = None
    usemotif = False
    if mode == 1:
        old_assign = np.loadtxt(oldAssignmentsName, dtype=int)
        usemotif = True
    (cluster_assignment, cluster_MRFs, motifs, motifRanked, bic, _) = solver.PerformFullTICC(
        initialClusteredPoints=old_assign, useMotif=usemotif)
    solver.CleanUp()
    if usemotif and outputDir is not None:
        # save the motifs and motifsRanked
        motifFile = "%smotifs.pkl" % outputDir
        pickleObject(motifFile, motifs)
        motifRankedName = "%smotifRanked.pkl" % outputDir
        pickleObject(motifRankedName, motifRanked)
    outputName = None
    if outputDir is not None:
        outputName = "%sassign.out" % outputDir
        np.savetxt(outputName, cluster_assignment, fmt='%d')
        outputName = "%sclusterMRFs" % outputDir
        pickleObject(outputName, cluster_MRFs)
    return outputName, bic

if __name__ == "__main__":
    # mode of 1 to skip old assign
    assert len(sys.argv) > 1
    mode = int(sys.argv[1])
    if mode == 2: 
        input_name = sys.argv[2]
        runBICTests(input_name, CLUSTER_NUMBER)
    else:
        assert len(sys.argv) == 4
        mode, input_fname, output_fdir = int(sys.argv[1]), sys.argv[2], sys.argv[3]
        dataset(mode, input_fname, output_fdir)


