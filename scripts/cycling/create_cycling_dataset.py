import csv
import numpy as np
import random
from collections import defaultdict
from sklearn.decomposition import PCA

# Dataset taken after downloading file directory from https://archive.ics.uci.edu/ml/datasets/Daily+and+Sports+Activities
PERSONS = ["1", "2", "3", "4", "5", "6", "7", "8"]
ACTIVITY = "15" # horizontal exercise bike
PERSON_MOTIF = ["1", "2", "3", "4"]


def read_block(filename):
    block = []
    with open(filename, "r") as f:
        csvreader = csv.reader(f, delimiter=",")
        for r in csvreader:
            block.append([float(x) for x in r])
    return block

def read_activity(activity, person):
    blocks = []
    for i in range(1, 61):
        data_name = str(i)
        if len(data_name) == 1:
            data_name = "0" + data_name
        activity_file_name = "data/a%s/p%s/s%s.txt" % (activity, person, data_name)
        blocks.append(read_block(activity_file_name))
    return blocks

activity_map = {p: read_activity(ACTIVITY, p) for p in PERSONS}

def pick_random_activity_block(person):
    block = random.choice(activity_map[person])
    labels = [person]*len(block)
    return labels, block

def pick_random_block():
    random_person = random.choice(PERSONS)
    return pick_random_activity_block(random_person)

NUM_GARBAGE = 4
def create_macro_segment():
    dataset = []
    labels = []
    motifMask = []
    for i in range(NUM_GARBAGE):
        person_labels, random_block = pick_random_block()
        labels += person_labels
        dataset += random_block
        motifMask += [0]*len(person_labels)
    for person in PERSON_MOTIF:
        person_labels, block = pick_random_activity_block(person)
        labels += person_labels 
        dataset += block
        motifMask += [1]*len(person_labels)
    return labels, dataset, motifMask

NUM_MACRO = 100
def create_entire_dataset():
    labels = []
    dataset = []
    motifMask = []
    for i in range(NUM_MACRO):
        seg_labels, seg_data, mm = create_macro_segment()
        labels += seg_labels
        dataset += seg_data
        motifMask += mm
    
    mean = np.mean(dataset, axis=0)
    dataset = dataset - mean
    stdev = np.std(dataset, axis=0)
    dataset = dataset/stdev 
    pca = PCA(n_components=10)
    dataset = pca.fit_transform(dataset)

    np.savetxt("activity_data.csv", dataset, delimiter=",")
    with open("labels.csv", "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(labels)
        # for l in labels:
        #     writer.writerow([l])
    with open("motifMask.csv", "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(motifMask)
    
    with open("labels_newlines.csv", "w") as csvfile:
        writer = csv.writer(csvfile)
        # writer.writerow(labels)
        for l in labels:
            writer.writerow([l])
    
    freqs = defaultdict(int)
    for w in labels:
        freqs[w] += 1
    print(freqs)

create_entire_dataset()

