import sys
import csv
from collections import Counter
import itertools
import heapq
import numpy as np
from sklearn.metrics import f1_score, confusion_matrix, accuracy_score
import pickle

NUM_CLUSTERS = 8


def performMapping(arr, perm):
    return [perm[i] for i in arr]


def getAssigns(fname):
    with open(fname, 'r') as instream:
        lines = instream.readlines()
        assigns = [int(val.strip()) for val in lines]
    return assigns

def getFullMapping(correctAssigns, testAssigns):
    ''' map test to correct '''
    mapping = getMappingDists(correctAssigns, testAssigns)
    # return performMapping(testAssigns, mapping)
    print(mapping)
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
    return mapping

def maskMotifSegments(correctAssigns, testAssigns, motifMask):
    correct = []
    test = []
    for i in range(len(motifMask)):
        if motifMask[i] != 0:
            correct.append(correctAssigns[i])
            test.append(testAssigns[i])
    return correct, test

def getValidMappings(correctAssigns, testAssigns, maskMotif):
    testMapped = getFullMapping(correctAssigns, testAssigns)
    # print(list(zip(correctAssigns, testMapped)))
    score = f1_score(correctAssigns, testMapped, average='weighted')
    acc = accuracy_score(correctAssigns, testMapped)
    print("total", score, acc)
    correctAssigns, testMapped = maskMotifSegments(correctAssigns, testMapped, maskMotif)
    # print(list(zip(correctAssigns, testMapped)))
    score = f1_score(correctAssigns, testMapped, average='weighted')
    acc = accuracy_score(correctAssigns, testMapped)
    print('motif', score, acc)

if __name__ == "__main__":
    assert len(sys.argv) == 4
    testfile = sys.argv[1]
    correct = sys.argv[2]
    maskMotif_file = sys.argv[3]

    testAssigns = getAssigns(testfile)
    print (len(testAssigns))
    with open(correct, "r") as correct:
        csvreader = csv.reader(correct)
        correctassigns = None
        for r in csvreader:
            correctassigns = [x for x in r]   
    correctmapping = list(set(correctassigns))
    correctmapping = {correctmapping[i]:i for i in range(len(correctmapping))}
    correctassigns = [correctmapping[c] for c in correctassigns]

    with open(maskMotif_file, "r") as f:
        csvreader = csv.reader(f)
        maskMotif = None
        for r in csvreader:
            maskMotif = [int(x) for x in r]    
    getValidMappings(correctassigns, testAssigns, maskMotif)

