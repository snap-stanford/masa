from CASC_solver import CASCSolver
import numpy as np
import sys
import pickle
import os
import time
CLUSTER_NUMBER = 10


def dataset(mode, input_Dir, motifReqs):
    beta = 25 # used 40 earlier
    number_of_clusters = CLUSTER_NUMBER
    oldAssignName = "%s/old/assign.out" % input_Dir
    input_name = "%s/data.out" % input_Dir
    if mode == 1: return runHyperParameterTests(input_name, input_Dir, number_of_clusters, beta, oldAssignName, motifReqs)
    return 0

def runHyperParameterTests(inputName, outputDir, clusters, beta, oldAssignmentsName, motifReqs):
    g = 0.8
    print(g, motifReqs)
    gammaDir = "%s/%s/" % (outputDir, g)
    makeDir(gammaDir)
    _, _, t = runTest(1, inputName, gammaDir, clusters,
                beta, g, motifReqs, oldAssignmentsName, 1)
    f = open("chkptTemp.out", "a+")
    f.write("%s, %s\n" % (gammaDir, t))
    return t

def runBICTests(inputName, number_of_clusters):
    beta = [25, 40, 60, 100]
    bicBeta = []
    for b in beta:
        _, bic = runNonMotifCASC(inputName, None, number_of_clusters, b, None)
        bicBeta.append((bic, b))
    bicBeta.sort(reverse=True)
    print(bicBeta)


def runNonMotifCASC(inputName, outputDir, clusters, beta, oldAssignmentsName):
    if outputDir is not None:
        oldDir = "%s/old/" % outputDir
        makeDir(oldDir)
        outputDir = oldDir
    return runTest(0, inputName, outputDir, clusters, beta, 1, 1, oldAssignmentsName, 15)

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
    print(oldAssignmentsName, inputName)
    solver = CASCSolver(window_size=1, number_of_clusters=clusters, lambda_parameter=1e-3, beta=beta, threshold=2e-5,
                        gamma=gamma, input_file=inputName, num_proc=30, maxMotifs=25, motifReq=motifReq, maxIters=maxIters)
    old_assign = None
    usemotif = False
    if mode == 1:
        old_assign = np.loadtxt(oldAssignmentsName, dtype=int)
        usemotif = True
    (cluster_assignment, cluster_MRFs, motifs, motifRanked, bic, t) = solver.PerformFullCASC(
        initialClusteredPoints=old_assign, useMotif=usemotif)
    solver.CleanUp()
    if usemotif and outputDir is not None:
        # save the motifs and motifsRanked
        motifFile = "%smotifs.pkl" % outputDir
        pickleObject(motifFile, motifs)
        motifRankedFile = "%smotifRanked.pkl" % outputDir
        pickleObject(motifRankedFile, motifRanked)
    outputName = None
    if outputDir is not None:
        outputName = "%sassign.out" % outputDir
        np.savetxt(outputName, cluster_assignment, fmt='%d')
    return outputName, bic, t

if __name__ == "__main__":
    #lenvals = ["250"]
    lenvals = ["1000", "1500", "2000", "3000", "4000", "5000", "6000", "7000"]
    tVals = []
    for lv in lenvals:
        directory = "scalability/%s/0.2" % lv
        tVals.append(dataset(1, directory, int(lv)/100))
    print(tVals, lenvals)
