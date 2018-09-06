import sys
import numpy as np
from sklearn.metrics import f1_score, accuracy_score
from generateDatasets.constants import NUM_CLUSTERS, GARBAGE_CLUSTERS, NUM_SEQS, NUM_GARBAGE, LEN_SEGMENT, CLUSTER_SEQUENCE
import pickle
import matplotlib.pyplot as plt
from collections import Counter


'''RUN  python3 analyze_synthetic_motifs.py ordered_synthetic_allperturbs/0.2/'''

GAMMAS = ["0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "0.99"]
ONLY_PICK_BEST = False

TOTAL_N = (NUM_GARBAGE + len(CLUSTER_SEQUENCE))*LEN_SEGMENT*NUM_SEQS

def getTrueResult():
    true_result = []
    for _ in range(NUM_SEQS):
        true_result += [0 for _ in range(NUM_GARBAGE*LEN_SEGMENT)]
        true_result += [1 for _ in range(len(CLUSTER_SEQUENCE)*LEN_SEGMENT)]
    return true_result

def getMotifResult(motifFname, motifRankedFname):
    f = open(motifFname, "rb")
    motifs = pickle.load(f) 

    f = open(motifRankedFname, "rb")
    motifRanked = pickle.load(f)

    motifResult = [0 for _ in range(TOTAL_N)]
    if ONLY_PICK_BEST:
        top = motifRanked[0][0]
        motifs  = {top:motifs[top]}
    for k,v in motifs.items():
        for s, e in v:
            motifResult[s:e+1] = [1 for _ in range(s, e+1)]
    return motifResult

def getBaselineScore(baseline="WEIGHTED"):
    true_result = getTrueResult()
    count = Counter(true_result)
    frac = float(count[1])/(count[0] + count[1])
    if baseline == "ALL":
        precision = frac
        recall = 1
        return (2*precision*recall/(precision+recall))
    elif baseline == "WEIGHTED":
        numOneCorrect = count[1]*frac # true positives
        numOneIncorrect = count[0]*frac # false positives
        precision = numOneCorrect/(numOneCorrect+numOneIncorrect)
        recall = numOneCorrect/(count[1])
        return ((2*precision*recall/(precision+recall)))
    elif baseline == "RANDOM":
        numOneCorrect = count[1]*0.5 # true positives
        numOneIncorrect = count[0]*0.5 # false positives
        precision = numOneCorrect/(numOneCorrect+numOneIncorrect)
        recall = numOneCorrect/(count[1])
        return ((2*precision*recall/(precision+recall)))
    else:
        assert False

def getScores(directory, baselineType="WEIGHTED"):
    cascscores = []
    trueResult = getTrueResult()
    baselineScore = getBaselineScore(baselineType)
    print(baselineType, baselineScore)
    baselines = [baselineScore for _ in range(len(GAMMAS))]
    labels = [float(v) for v in GAMMAS]
    for g in GAMMAS:
        motifName = "%s/%s/motifs.pkl" % (directory, g)
        motifRanksName = "%s/%s/motifRanked.pkl" % (directory, g)
        motifResult = getMotifResult(motifName, motifRanksName)
        cascscores.append(f1_score(trueResult, motifResult))

    plt.figure(1, figsize=(7, 4))
    print(cascscores)
    plt.plot(GAMMAS, cascscores, '-bv', label='CASC')
    # plt.plot(GAMMAS, baselines, '--', label='Random', color="orange")
    # plt.plot(GAMMAS, all_scores[2], '-.', label='No motif')
    plt.xlabel("$\gamma$")
    plt.ylim(ymin=0.2, ymax=1)
    plt.ylabel("F1 Score")
    # plt.title("Motif Identification F1 vs $\gamma$, $\epsilon = 0.4$")
    plt.legend(loc='lower center', ncol=2, fancybox=False, edgecolor="black")
    plt.show()


    # np.savetxt(outputfile, all_scores, delimiter=",")

if __name__ == "__main__":
    fileDir = sys.argv[1]
    getScores(fileDir, "WEIGHTED")
