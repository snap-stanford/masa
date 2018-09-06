import numpy as np
from snap import GenRndGnm, PNGraph

SAVE_COVS = False


def genInvCov(size, low=0.3, upper=0.6, portion=0.2):
    ''' Generate off diagonal blocks'''
    portion = portion/2
    S = np.zeros((size, size))
    G = GenRndGnm(PNGraph, size, int((size*(size-1))*portion))
    for EI in G.Edges():
        value = (np.random.randint(2) - 0.5)*2 * \
            (low + (upper - low)*np.random.rand(1)[0])
        S[EI.GetSrcNId(), EI.GetDstNId()] = value
    S = S + S.T
    return np.matrix(S)


def genRandInv(size, low=0.3, upper=0.6, portion=0.2):
    ''' Generate diagonal blocks '''
    S = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            if np.random.rand() < portion:
                value = (np.random.randint(2) - 0.5)*2 * \
                    (low + (upper - low)*np.random.rand(1)[0])
                S[i, j] = value
    return np.matrix(S)


def generate_inverse(n, sparsity, w):
    '''
    n: number of sensors
    w: window size
    sparsity: sparsity penalty (lambda)
    '''
    N = n*w
    block_matrices = {}
    # Generate all the blocks
    # create diagonal blocks
    block_matrices[0] = genInvCov(size=n, portion=sparsity)
    for block in range(1, w):
        block_matrices[block] = genRandInv(size=n, portion=sparsity)

    inv_matrix = np.zeros([N, N])
    # go through all the blocks
    for i in range(w):
        for j in range(w):
            block_result = block_matrices[np.abs(i - j)]
            if i <= j:
                block_result = block_result.T
            inv_matrix[i*n:(i+1)*n, j*n:(j+1)*n] = block_result

    # Make the matrix positive definite
    lambda_min = np.min(np.linalg.eig(inv_matrix)[0])

    inv_matrix = inv_matrix + (0.1 + abs(lambda_min)) * np.identity(N)
    eigs, _ = np.linalg.eig(inv_matrix)
    return inv_matrix


def getMeanCov(dataCounter, cluster, cluster_covs, cluster_mean, cluster_mean_stacked, Data, w, n, K):
    N = n*w
    if dataCounter == 0:
        cov_matrix = cluster_covs[cluster][0:n, 0:n]
        new_mean = cluster_mean_stacked[n * (w-1):N]
        return new_mean.reshape(n), cov_matrix
    breakIndex = dataCounter  # update dataCounter
    # generate the following points
    num = breakIndex if breakIndex < w else w-1
    cov_matrix = cluster_covs[cluster][0:(num+1)*n, 0:(num+1)*n]
    Sig22 = cov_matrix[(num)*n:(num+1)*n, (num)*n:(num+1)*n]
    Sig11 = cov_matrix[0:(num)*n, 0:(num)*n]
    Sig21 = cov_matrix[(num)*n:(num+1)*n, 0:(num)*n]
    Sig12 = Sig21.T
    # sigma2|1
    sig21_invsig11 = np.dot(Sig21, np.linalg.inv(Sig11))
    cov_mat_tom = Sig22 - np.dot(sig21_invsig11, Sig12)
    a = np.zeros((num*n, 1))
    for idx in range(num):
        if breakIndex < w:
            a[idx*n:(idx+1)*n, 0] = Data[idx, :].reshape((n))
        else:
            a[idx*n:(idx+1)*n, 0] = Data[breakIndex -
                                         w + 1 + idx, :].reshape([n])
    new_mean = cluster_mean + \
        np.dot(sig21_invsig11, (a - cluster_mean_stacked[0:(num)*n, :]))
    return new_mean.reshape(n), cov_mat_tom


def generate_data(K, n, w, sparsity, assignments, out_file_name, noiseWeight=1):
    '''
    K: number of clusters
    n: number of sensors
    w: window size
    assignments: [(cluseterID1, clusterID2), numberOfPoints)]
        if clusterID2 = None, then we're not adding noise
    out_file_name: outputfile to save assignments
    noiseWeight: how much weight should clusterID2 get
    '''

    # GENERATE POINTS
    N = n*w
    cluster_mean = np.zeros([n, 1])
    cluster_mean_stacked = np.zeros([N, 1])

    # Generate two inverse matrices
    cluster_invs = {}
    cluster_covs = {}
    for c in range(K):
        cluster_invs[c] = generate_inverse(n, sparsity, w)
        cluster_covs[c] = np.linalg.inv(cluster_invs[c])
        if SAVE_COVS:
            fn = "inv_cov_cluster_%s.csv" % (c)
            np.savetxt(fn, cluster_invs[c], delimiter=",", fmt='%1.6f')
            fn = "cov_cluster_%s.csv" % (c)
            np.savetxt(fn, cluster_covs[c], delimiter=",", fmt='%1.6f')

    # Data matrix
    T = sum([val[1] for val in assignments])
    print(T)
    Data = np.zeros((T, n))
    dataCounter = 0
    for clusterList, numPoints in assignments:
        cluster1, cluster2 = clusterList
        for i in range(numPoints):
            mean, cov = None, None
            mean1, cov1 = getMeanCov(dataCounter, cluster1, cluster_covs,
                       cluster_mean, cluster_mean_stacked, Data, w, n, K)
            if cluster2 is None:
                mean, cov = mean1, cov1 
            else:
                mean2, cov2 = getMeanCov(dataCounter, cluster2, cluster_covs,
                       cluster_mean, cluster_mean_stacked, Data, w, n, K)
                mean = (1-noiseWeight)*mean1 + noiseWeight*mean2
                cov = (1-noiseWeight)*cov1 + noiseWeight*cov2
            new_row = np.random.multivariate_normal(mean, cov)
            Data[dataCounter, :] = new_row
            dataCounter += 1

    print("length of generated Data is:", Data.shape[0])
    np.savetxt(out_file_name, Data, delimiter=",", fmt='%1.4f')
