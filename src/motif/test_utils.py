from hmm import MotifHMM
import numpy as np

def GenerateFakeData(assignments, confidenceMean, confidenceVariance, numClusters):
    '''
    Create a probability matrix given the above breakpoints
    @ params
    breakpoints - a list of the length of each cluster, with the first and last values being random
    '''
    likelihoods = []
    for c in assignments:
        if c == -1:
            # this is garbage
            a = np.random.rand(numClusters)
            a = a/np.sum(a)
            likelihoods.append(a)
        else:
            a = np.random.rand(numClusters)
            a[c] = 0
            val = min(0.99, np.random.normal(confidenceMean, confidenceVariance))
            rest = 1-val
            a = rest*a/np.sum(a)
            a[c] = val
            likelihoods.append(a)
    likelihoods = -1*np.log(np.array(likelihoods))
    return likelihoods

def test(dataAssign, motif, confidenceMean, confidenceVariance, numClusters):
    ll = GenerateFakeData(dataAssign, confidenceMean, confidenceVariance, numClusters)
    motif_hmm = MotifHMM(ll, motif,0.9, None)
    result,_ = motif_hmm.SolveAndReturn()
    ll = motif_hmm.negLLMatrix
    for i in range(len(dataAssign)):
        print "%s \t %s \t %s" % (dataAssign[i], result[i], ll[i])


def testAssign1():
    clusterLen = 5
    dataAssign = [-1]*clusterLen + [0]*clusterLen + [1]*clusterLen + [2]*clusterLen + [-1]*clusterLen  + [-1]*clusterLen + [0]*clusterLen + [1]*clusterLen + [2]*clusterLen + [-1]*clusterLen
    motif = [0, 1, 2]
    confidenceMean = 0.5
    confidenceVariance = 0.05
    numClusters = 4
    test(dataAssign, motif, confidenceMean, confidenceVariance, numClusters)

testAssign1()