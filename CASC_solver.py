import numpy as np
import collections
import logging

from scipy import stats

from sklearn import mixture
import time
from multiprocessing import Pool
from collections import deque
from src.CASC_helper import *
from src.motif.find_motif import PerformAssignment
from src.admm_solver import ADMMSolver
#######################################################################################################################################################################
np.set_printoptions(formatter={'float': lambda x: "{0:0.4f}".format(x)})
np.random.seed(103)  # 102

#####################################################################################################################################################################################################

logging.basicConfig(level=logging.DEBUG)


class CASCSolver:
    def __init__(self, window_size=10, number_of_clusters=5, lambda_parameter=11e-2,
                 beta=400, maxIters=1000, threshold=2e-5,
                 input_file=None, num_proc=1, gamma=0.9, maxMotifs=None, motifReq=2):
        self.window_size = window_size
        self.K = number_of_clusters  # number of clusters
        self.lambda_param = lambda_parameter
        self.beta = beta
        self.maxIters = maxIters
        self.threshold = threshold
        self.num_proc = num_proc
        self.cluster_reassignment = 20  # number of points to reassign to a 0 cluster
        self.gamma = gamma  # aggressiveness of motif
        self.maxMotifs = maxMotifs  # max num of motifs
        self.motifReq = motifReq
        # get the data inflated by window size
        data = np.loadtxt(input_file, delimiter=",")
        m, n = data.shape
        self.m = m  # observations
        self.n = n  # size of each observation vector
        logging.info("done retrieving data")
        self.complete_data = np.zeros([m, window_size*n])
        for i in range(m):
            for k in range(window_size):
                if i+k < m:
                    self.complete_data[i][k*n:(k+1)*n] = data[i+k][0:n]

        # start the process pool
        self.pool = Pool(processes=self.num_proc)

    def CleanUp(self):
        if self.pool is not None:
            self.pool.close()
            self.pool.join()

    def PerformFullCASC(self, initialClusteredPoints=None, useMotif=False):
        train_cluster_inverse = motifs = motifsRanked = None
        clustered_points = initialClusteredPoints
        start = time.time()
        if clustered_points is None:
            # perform no motif CASC
            initialClusteredPoints = self.getInitialClusteredPoints()
            clustered_points, train_cluster_inverse, _, _, bic = self.solveWithInitialization(
                initialClusteredPoints, useMotif=False)
        if useMotif:
            # perform secondary motif CASC if specified
            clustered_points, train_cluster_inverse, motifs, motifsRanked, bic = self.solveWithInitialization(
                clustered_points, useMotif=True)
        end = time.time()
        return clustered_points, train_cluster_inverse, motifs, motifsRanked, bic, end - start

    def getInitialClusteredPoints(self):
        gmm = mixture.GaussianMixture(
            n_components=self.K, covariance_type="full")
        gmm.fit(self.complete_data)
        clustered_points = gmm.predict(self.complete_data)
        return clustered_points

    def solveWithInitialization(self, clustered_points, useMotif):
        motifs = None
        rankedMotifs = None
        assert self.maxIters > 0  # must have at least one iteration
        num_stacked = self.window_size
        beta = self.beta  # switching penalty
        lam_sparse = self.lambda_param  # sparsity parameter
        K = self.K  # Number of clusters

        logging.info("lambda: %s, beta: %s, clusters: %s, num stacked %s" % (
            lam_sparse, beta, K, num_stacked))

        train_cluster_inverse = {}
        computed_cov = {}
        cluster_mean_stacked_info = {}
        clustered_point_historyN = 3
        clustered_point_history = deque(
            [None for i in range(clustered_point_historyN)])
        empirical_covariances = {}

        # PERFORM TRAINING ITERATIONS
        for iters in range(self.maxIters):
            logging.debug("\n\n\n ITERATION ### %s" % iters)

            # {cluster: [point indices]}
            clust_indices = self.getClustIndices(clustered_points)

            # solve for clusters
            self.solveForClusters(clust_indices, cluster_mean_stacked_info,
                                  empirical_covariances, train_cluster_inverse, computed_cov)

            # update old computed covariance
            old_computed_cov = computed_cov

            LLE_all_points_clusters = self.getLikelihood(
                computed_cov, cluster_mean_stacked_info, clustered_points)

            # Update cluster points
            clustered_points = updateClusters(
                LLE_all_points_clusters, switch_penalty=beta)
            clustered_points = [i.astype(int) for i in clustered_points]
            if useMotif:
                clustered_points, motifs, rankedMotifs = PerformAssignment(
                    clustered_points, LLE_all_points_clusters, self)
                for m, score in rankedMotifs:
                    print("%s ---> %s, %s" % (m, len(motifs[m]), score))
            before_zero = clustered_points.copy()
            self.assignToZeroClusters(
                clustered_points, old_computed_cov, computed_cov, cluster_mean_stacked_info)
            logging.debug("smoothened points")

            clust_indices = self.getClustIndices(clustered_points)
            for cluster in range(K):
                logging.debug(
                    "length of cluster %s --> %s" % (cluster, len(clust_indices[cluster])))
            stop_here = False
            if before_zero in clustered_point_history:
                logging.info("CONVERGED!!!! BREAKING EARLY!!!")
                stop_here = True
            clustered_point_history.popleft()
            clustered_point_history.append(before_zero)
            if stop_here:
                break
        bic = computeBIC(self.K, self.m, clustered_point_history[-1], train_cluster_inverse,
                         empirical_covariances)
        clusterbic = None
        if motifs is not None:
            clusterbic = computeClusterBIC(
                self.K, self.m, clustered_points, train_cluster_inverse, empirical_covariances, motifs)
        logging.info("BIC for beta %s clusters %s is %s" % (self.beta, self.K, bic))
        return (clustered_point_history[-1], train_cluster_inverse, motifs, rankedMotifs, (bic, clusterbic))

    def getLikelihood(self, computed_cov, cluster_mean_stacked_info, clustered_points):
        '''
        Get the likelihood matrix
        '''
        K = self.K
        num_blocks = self.window_size + 1
        num_stacked = self.window_size
        n = self.n
        N = len(clustered_points)
        inv_cov_dict = {}  # cluster to inv_cov
        log_det_dict = {}  # cluster to log_det
        for cluster in range(K):
            cov_matrix = computed_cov[cluster][0:(
                num_blocks-1)*n, 0:(num_blocks-1)*n]
            inv_cov_dict[cluster] = np.linalg.inv(cov_matrix)
            log_det_dict[cluster] = np.log(np.linalg.det(cov_matrix))

        LLE_all_points_clusters = np.zeros([N, K])
        for point in range(N):
            if point + num_stacked-1 < self.complete_data.shape[0]:
                for cluster in range(K):
                    cluster_mean_stacked = cluster_mean_stacked_info[cluster]
                    x = self.complete_data[point, :] - \
                        cluster_mean_stacked[0:(num_blocks-1)*n]
                    inv_cov_matrix = inv_cov_dict[cluster]
                    log_det_cov = log_det_dict[cluster]
                    lle = np.dot(x.reshape([1, (num_blocks-1)*n]), np.dot(
                        inv_cov_matrix, x.reshape([n*(num_blocks-1), 1]))) + log_det_cov
                    LLE_all_points_clusters[point, cluster] = lle
        # normalize
        ll = -1*LLE_all_points_clusters
        '''
        normalizer = np.reshape(np.max(ll, axis=1), (ll.shape[0], 1))
        ll = ll - normalizer
        '''
        ll = np.exp(ll)
        '''
        sums = np.reshape(np.sum(ll, axis=1),(ll.shape[0], 1))
        ll = ll/sums
        '''
        ll[ll < 1e-320] = 1e-320
        assert np.all(ll > 0)
        ll = -1*np.log(ll)
        return ll

    def assignToZeroClusters(self, clustered_points, old_computed_cov, computed_cov, cluster_mean_stacked_info):
        '''
        N should be length of clustered_points
        '''
        K = self.K
        clust_indices = self.getClustIndices(clustered_points)
        N = len(clustered_points)
        cluster_lens = {k: len(clust_indices[k]) for k in range(K)}
        cluster_norms = [(np.linalg.norm(old_computed_cov[i]), i)
                         for i in range(K)]
        norms_sorted = sorted(cluster_norms, reverse=True)
        # clusters that are not 0 as sorted by norm
        valid_clusters = [cp[1]
                          for cp in norms_sorted if cluster_lens[cp[1]] != 0]

        # Add a point to the empty clusters
        # assuming more non empty clusters than empty ones
        counter = 0
        for cluster in range(K):
            if cluster_lens[cluster] == 0:
                # a cluster that is not len 0
                cluster_selected = valid_clusters[counter]
                counter = (counter+1) % len(valid_clusters)
                # random point number from that cluster
                start_point = np.random.choice(clust_indices[cluster_selected])
                for i in range(0, self.cluster_reassignment):
                    # put cluster_reassignment points from point_num in this cluster
                    point_to_move = start_point + i
                    if point_to_move >= N:
                        break
                    # update stats
                    clustered_points[point_to_move] = cluster
                    computed_cov[cluster] = old_computed_cov[cluster_selected]
                    cluster_mean_stacked_info[cluster] = self.complete_data[point_to_move, :]

    def solveForClusters(self, clust_indices, cluster_mean_stacked_info,
                         empirical_covariances, train_cluster_inverse, computed_cov):
        '''
        Find the characteristic clusters. Given clust_indices, fill out the results
        in the rest of the parameters
        '''
        K = self.K
        num_stacked = self.window_size
        n = self.n
        cluster_lens = {k: len(clust_indices[k]) for k in range(K)}

        optRes = [None] * K
        for cluster in range(K):
            if cluster_lens[cluster] != 0:
                cluster_data = np.take(
                    self.complete_data, clust_indices[cluster], axis=0)
                cluster_mean_stacked_info[cluster] = np.mean(
                    cluster_data, axis=0)
                # Fit a model - OPTIMIZATION
                probSize = num_stacked * n
                lamb = np.zeros((probSize, probSize)) + self.lambda_param
                S = np.cov(np.transpose(cluster_data))
                empirical_covariances[cluster] = S
                solver = ADMMSolver(lamb, num_stacked, n, 1, S)
                # apply to process pool
                optRes[cluster] = self.pool.apply_async(
                    solver, (1000, 1e-6, 1e-6, False,))

        for cluster in range(K):
            if optRes[cluster] == None:
                continue
            val = optRes[cluster].get()
            logging.debug("optimization for cluster %s is done" % cluster)
            S_est = upperToFull(val, 0)
            X2 = S_est
            cov_out = np.linalg.inv(X2)
            computed_cov[cluster] = cov_out
            train_cluster_inverse[cluster] = X2

    def getClustIndices(self, clustered_points):
        clust_indices = collections.defaultdict(list)
        for point, cluster in enumerate(clustered_points):
            clust_indices[cluster].append(point)
        return clust_indices
