from generate_synthetic_util import generate_data
import numpy as np
import sys
from constants import NUM_CLUSTERS, GARBAGE_CLUSTERS, CLUSTER_SEQUENCE, NUM_SEQS, NUM_GARBAGE, LEN_SEGMENT

DELTA = 0.6
WINDOW_SIZE = 1
NUM_SENSORS = 5
SPARSITY = 0.2
RAND_SEED = 30 # 30
np.random.seed(RAND_SEED)

'''
NUM_CLUSTERS = 10
GARBAGE_CLUSTERS = 10
CLUSTER_SEQUENCE = [6,7,8,9]
NUM_SEQS = 500 # number of macro segs
NUM_GARBAGE = 10# number of garbage segs
LEN_SEGMENT = 25 # length of each segment
'''

def createCorrect(outputFilename):
    ''' just creates in terms of cluster segments'''
    assigns = []
    for _ in range(NUM_SEQS):
        assigns += np.random.choice(GARBAGE_CLUSTERS, NUM_GARBAGE).tolist()
        assigns += CLUSTER_SEQUENCE
    # pre-pend with length
    np.savetxt(outputFilename, np.array(assigns), delimiter=",", fmt='%d')

def createDataset(correctFileName, eps, outputFilename):
    # delta: weight of other cluster that is being used to perturb
    # eps: fraction of motif segments being perterubed
    # rho: integer <= 4 which is number of legs we are perturbing
    # correctFilename: the ground truth
    delta = DELTA

    groundTruth = []
    with open(correctFileName, 'r') as instream:
        lines = instream.readlines()
        groundTruth = [int(val.strip()) for val in lines]

    assigns = []
    lenMacro = len(CLUSTER_SEQUENCE) + NUM_GARBAGE
    for i in range(NUM_SEQS):
        segment = []
        macroSegment = groundTruth[i*lenMacro: (i+1)*lenMacro]
        for clust in macroSegment:
            otherClust = None
            if np.random.random() < eps:
                otherClust = np.random.choice(CLUSTER_SEQUENCE[0])
            segment += [(clust, otherClust)]
        assigns += segment

#         garbageSegment = macroSegment[:NUM_GARBAGE]
#         perturbs = macroSegment[NUM_GARBAGE:]
#         segment += [(a, None) for a in garbageSegment]
#         if np.random.random() > eps or rho == 0:
#             segment += [(a, None) for a in perturbs]
#         else:
#             perturbIndices = np.random.choice(np.arange(len(CLUSTER_SEQUENCE)), rho, replace=False)
#             for p in range(len(CLUSTER_SEQUENCE)):
#                 if p in perturbIndices:
#                     # perturb this leg
#                     #l = [num for num in range(NUM_CLUSTERS) if num != perturbs[p]]
#                     chosenCluster = np.random.choice(CLUSTER_SEQUENCE[0])
#                     segment += [(perturbs[p], chosenCluster)]
#                 else: segment += [(perturbs[p], None)]
#         #print(segment)
#         assigns += segment
    assignment, _ = createSegments(assigns, LEN_SEGMENT)
    print(delta)
    generate_data(NUM_CLUSTERS, NUM_SENSORS, WINDOW_SIZE,
                  SPARSITY, assignment, outputFilename, noiseWeight=delta)


def createSegments(assigns, lenSegment):
    assignment = []
    correctAssignment = []
    for cluster, othercluster in assigns:
        correctAssignment += np.full(lenSegment, cluster, dtype=int).tolist()
        assignment.append(((cluster, othercluster), lenSegment))
    return assignment, correctAssignment


if __name__ == "__main__":
    assert len(sys.argv) > 2
    mode = int(sys.argv[1])
    if mode == 0:
        outputFile = sys.argv[2]
        createCorrect(outputFile)
    else:
        assert len(sys.argv) == 4
        # correctFileName, delta, eps, rho, outputFilename
        outputDir, e = sys.argv[2], float(sys.argv[3])
        corrfname = "%s/correct.out" % outputDir
        outputFilename = "%s/%s/data.out" % (outputDir, sys.argv[3])
        print(outputFilename)
        createDataset(corrfname, e, outputFilename)
