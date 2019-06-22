from TICC_solver import TICCSolver
import numpy as np
import sys

assert len(sys.argv) == 7


def test():
    mode = int(sys.argv[1])
    clusters = int(sys.argv[2])
    beta = float(sys.argv[3])
    inputName = sys.argv[4]
    old_assignmentsName = sys.argv[5]
    outputName = sys.argv[6]
    if mode == 1:
        runHyperParameterTests(inputName, outputName,
                               clusters, beta, old_assignmentsName)
    else:
        runNonMotifTICC(inputName, outputName, clusters,
                        beta, old_assignmentsName)


def runHyperParameterTests(inputName, outputName, clusters, beta, oldAssignmentsName):
    gammas = [0.6, 0.7, 0.8, 0.9, 0.99]
    motifReqs = [3]
    motifDict = {}
    for g in gammas:
        for m in motifReqs:
            motifs, rankedMotifs = runTest(1, inputName, outputName, clusters,
                             beta, g, m, oldAssignmentsName)
            motifDict[(g, m)] = (motifs, rankedMotifs)
    print(motifDict)


def runNonMotifTICC(inputName, outputName, clusters, beta, oldAssignmentsName):
    runTest(0, inputName, outputName, clusters, beta, 1, 1, oldAssignmentsName)


def runTest(mode, inputName, outputName, clusters, beta, gamma, motifReq, oldAssignmentsName):
    print("TESTING %s %s" % (beta, gamma))
    solver = TICCSolver(window_size=3, number_of_clusters=clusters, lambda_parameter=5e-3, beta=beta, threshold=2e-5,
                        gamma=gamma, input_file=inputName, num_proc=30, maxMotifs=50, motifReq=motifReq, maxIters=50)
    old_assign = None
    usemotif = False
    if mode == 1:
        print("using motif")
        old_assign = np.loadtxt(oldAssignmentsName, dtype=int)
        usemotif = True

    (cluster_assignment, cluster_MRFs, motifs, rankedMotifs) = solver.PerformFullTICC(
        initialClusteredPoints=old_assign, useMotif=usemotif)
    solver.CleanUp()
    if mode == 1:
        fname = "%s_clust%s_beta%s_gamma%s_req%s.out" % (
            outputName, clusters, beta, gamma, motifReq)
    else:
        fname = "%s_clust%s_beta%s.out" % (outputName, clusters, beta)
    print ("saving to ", fname)
    np.savetxt(fname, cluster_assignment, fmt='%d')
    return motifs, rankedMotifs


test()
