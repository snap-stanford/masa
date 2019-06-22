import numpy as np 
from sklearn.cluster import KMeans
from hmmlearn.hmm import GaussianHMM
from sklearn.mixture import GaussianMixture
import sys

def getData(data_file):
    data = np.loadtxt(data_file, delimiter=",")
    return data

def performBaseline(data_file, out_file, baselineType="KMEANS", K=10):
    data = np.loadtxt(data_file, delimiter=",")
    labels=None
    if baselineType == "KMEANS":
        kmeans = KMeans(n_clusters=K, random_state=0, n_jobs=30)
        kmeans.fit(data)
        labels = kmeans.labels_
    elif baselineType == "GMM":
        gmm = GaussianMixture(n_components=K, covariance_type="full")
        gmm.fit(data)
        labels = gmm.predict(data)
    elif baselineType == "HMM":
        hmm = GaussianHMM(n_components=K, n_iter=100, random_state=100)
        hmm.fit(data)
        labels = hmm.predict(data)
    else:
        assert False

    labels = np.array(labels)
    np.savetxt(out_file, labels, delimiter=",", fmt='%d')

if __name__ == "__main__":
    assert len(sys.argv) > 1
    directory = sys.argv[1]
    #mapping = {"HMM": "hmm.out"}
    mapping = {"KMEANS":"kmeans.out", "GMM": "gmm.out", "HMM": "hmm.out"}
    #epsilons = ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6"]
    epsilons = ["0.05"]
    #epsilons = ["0.0", "0.1", "0.2"]
    for e in epsilons:
        infile = "%s/%s/data.out" % (directory, e)
        for k,v in mapping.items():
            outname = "%s/%s/%s" % (directory, e, v)
            print(infile, outname, k)
            performBaseline(infile, outname, k, 10)
    '''
    assert len(sys.argv) > 3
    baseline = sys.argv[1]
    infile = sys.argv[2]
    outfile = sys.argv[3]
    performBaseline(infile, outfile, baseline)
    '''
