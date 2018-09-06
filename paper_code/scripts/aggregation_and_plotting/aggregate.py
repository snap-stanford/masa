import sys
"""Collapse a CASC output"""
fname = sys.argv[1]

with open(fname, "rb") as file:
    currIndex = -1
    index = -1
    currCluster = -1
    for line in file:
        cluster = int(line)
        if currCluster != cluster:
            if currCluster != -1:
                print currIndex, index, currCluster
            index += 1
            currIndex = index
            currCluster = cluster
        else:
            index += 1
    print "%d,%d,%d" % (currIndex, index, cluster)
