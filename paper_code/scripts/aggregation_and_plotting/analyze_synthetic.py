import sys
from collections import Counter
import itertools
import heapq
import numpy as np
from sklearn.metrics import f1_score, confusion_matrix, accuracy_score
import matplotlib
matplotlib.use('agg')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from generateDatasets.constants import NUM_CLUSTERS, GARBAGE_CLUSTERS, NUM_SEQS, NUM_GARBAGE, LEN_SEGMENT, CLUSTER_SEQUENCE
rcParams.update({'figure.autolayout': True})
import pickle

FIGURE_COUNT = 1
PLOT_CONF = False
VERBOSE = False
STARTING_MOTIF = CLUSTER_SEQUENCE[0] # i.e 6


def plot_confusion_matrix(cm, classes, output_name,
                          normalize=False,
                          title='Confusion matrix',
                          cmap=plt.cm.Blues,):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    # print(cm)

    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt), fontsize=10,
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    # plt.show()


def performMapping(arr, perm):
    return [perm[i] for i in arr]

def plotTable(text):
    # global FIGURE_COUNT
    # plt.figure(FIGURE_COUNT)
    # FIGURE_COUNT += 1
    fig, ax = plt.subplots()
    # Hide axes
    ax.xaxis.set_visible(False) 
    ax.yaxis.set_visible(False)

    # Table from Ed Smith answer
    collabel=list(range(NUM_CLUSTERS))
    rowlabel = list(range(NUM_CLUSTERS))
    ax.table(cellText=text,colLabels=collabel,loc='center', rowLabels=rowlabel)
    plt.show()

def getAssigns(fname):
    with open(fname, 'r') as instream:
        lines = instream.readlines()
        assigns = [int(val.strip()) for val in lines]
    return assigns

def getMask(assign):
    lenBlock = (NUM_GARBAGE + len(CLUSTER_SEQUENCE))*LEN_SEGMENT
    lenGarbage = NUM_GARBAGE*LEN_SEGMENT
    numSeqs = NUM_SEQS
    result = []
    for s in range(numSeqs):
        result += assign[s*lenBlock + lenGarbage :s*lenBlock + lenBlock + 1]
    return result

def getFullMapping(correctAssigns, testAssigns):
    ''' map test to correct '''
    assignCorrectMask = getMask(correctAssigns[:])
    assignTestMask = getMask(testAssigns[:])
    mapping = getMappingDists(assignCorrectMask, assignTestMask)
    mapping[:STARTING_MOTIF] = [None for _ in range(STARTING_MOTIF)]
    mapping = getMappingDists(correctAssigns, testAssigns, mapping)
    reverseMapping = [None for _ in range(len(mapping))]
    for i,val in enumerate(mapping):
        reverseMapping[val] = i
    return performMapping(testAssigns, reverseMapping)

def getMappingDists(correctAssigns, testAssigns, mapping=None):
    ''' mapping from testAssigns to correctAssigns'''
    T = len(correctAssigns)
    K = NUM_CLUSTERS
    results = [np.zeros(K) for i in range(K)]
    for i in range(T):
        correct = correctAssigns[i]
        test = testAssigns[i]
        results[correct][test] += 1.0
    x = []
    indexes = [i for i in range(K)]
    #plotTable(results)
    for i in range(K):
        final_result = results[i]
        final_result *= -1
        # final_result /= -1*np.sum(results[i]) # invert order
        final_result = final_result.tolist()
        indices_result = list(zip(final_result, indexes[:]))
        indices_result.sort()
        x.append((indices_result, i))
        if VERBOSE:
            resultStr = ["%s:%d" % (j, results[i][j]) for j in range(K)]
            print("%s: %s" % (i, resultStr))
    heapq.heapify(x)
    taken = [False for _ in range(K)]
    if mapping is None:
        mapping = [None for _ in range(K)]
    for val in mapping:    
        if val is not None:
           taken[val] = True
    while len(x) != 0:
        r, idx = heapq.heappop(x)
        if mapping[idx] is not None: continue
        best_score, best_index = r[0]
        if taken[best_index]:
            r = r[1:]
            assert len(r) != 0
            heapq.heappush(x, (r, idx))
        else:
            mapping[idx] = best_index
            taken[best_index] = True
    if VERBOSE: print(mapping)
    return mapping


def getValidMappings(correctFname, assignFname):
    correctAssignIndices = getAssigns(correctFname)
    correctAssigns = []
    for a in correctAssignIndices: correctAssigns += [a for _ in range(LEN_SEGMENT)]
    testAssigns = getAssigns(assignFname)
    testMapped = getFullMapping(correctAssigns, testAssigns)
    cared_about_values = list(range(STARTING_MOTIF,NUM_CLUSTERS))
    score = f1_score(correctAssigns, testMapped,
                     average='weighted', labels=cared_about_values)
    score2 = f1_score(correctAssigns, testMapped, average='weighted')
    print(cared_about_values)
    correctMotifMask = getMask(correctAssigns)
    testMotifMask = getMask(testMapped)
    score3 = accuracy_score(correctMotifMask, testMotifMask)
    cf = confusion_matrix(correctAssigns, testMapped)
    global FIGURE_COUNT, PLOT_CONF
    if PLOT_CONF:
        plt.figure(FIGURE_COUNT)
        FIGURE_COUNT += 1
        plot_confusion_matrix(cf, list(range(NUM_CLUSTERS)),"blah", title='Confusion matrix specific', normalize=True)
    #print("---")
    print("only relevant: %s, all: %s, accuracy: %s" % (score, score2, score3))
    return score, score2, score3

if __name__ == "__main__":
    assert len(sys.argv) == 3
    testfile = sys.argv[1]
    correct = sys.argv[2]
    getValidMappings(correct, testfile)
if PLOT_CONF:
    plt.show()

