import numpy as np

# TODO Add duration modelling (honestly v necesary)


class HMM:
    '''
    '''

    def __init__(self, adjacencyMatrix, negLLMatrix, initDistribution):
        '''
        @params for m = # states, n = # data points
        adjacencymatrix: mxm matrix where A[i,j] is the switching penalty of going 
            from state i to state j (infty if not possible). 
        negLLMatrix: nxm matrix where A[i, j] is the negative log likelihood of 
            the time data point at i being in state j
        initDistribution: mx1 matrix the neg ll distribution of initial states
            (i.e, infty if you don't want a state to be first)
        '''
        self.n, self.m = np.shape(negLLMatrix)
        self.adjacencyMatrix = adjacencyMatrix
        self.negLLMatrix = negLLMatrix
        self.initDistribution = initDistribution
        self.viterbiGrid = np.zeros((self.n, self.m))  # dp grid
        self.backPointers = np.full((self.n, self.m), -1)  # backpointer grid

    def GetAdjacencyMatrix(self, ts):
        ''' overridden by motif hmm'''
        return self.adjacencyMatrix

    def UpdateStep(self, ts):
        '''
        Update the ts column of the viterbi grid
        '''
        if ts == 0:  # this is the first column
            self.viterbiGrid[0] = self.initDistribution + self.negLLMatrix[0]
            return
        # the previous costs
        prevCosts = np.reshape(self.viterbiGrid[ts - 1], (self.m, 1))
        # compute costs going to j assuming you are starting at i
        withSwitchingCosts = prevCosts + self.GetAdjacencyMatrix(ts)
        indices = np.argmin(withSwitchingCosts, axis=0)  # backpointers
        maxvals = withSwitchingCosts[indices, range(self.m)]
        indices[maxvals == np.infty] = -1
        newRow = maxvals + self.negLLMatrix[ts]  # the new costs
        self.viterbiGrid[ts] = newRow
        self.backPointers[ts] = indices

    def GenerateSequenceFromBackPointer(self, ending=None):
        '''
        Starting from the end to 0, construct the most likely sequence.
        When ending is not null, force the last stage to be ending
        '''
        # ts is the ending index
        seq = [ending]
        currStage = ending
        endingIndex = self.n - 1
        for i in range(endingIndex, 0, -1):
            # going backward, skipping first pointer
            currStage = int(self.backPointers[i][currStage])
            seq = [currStage] + seq
        return seq

    def solveHMM(self):
        for i in range(self.n):
            self.UpdateStep(i)

    def getEndingScore(self, state):
        ''' get the last score for the given state'''
        return self.viterbiGrid[-1, state]


class MotifHMM(HMM):
    def __init__(self, negLLMatrix, motif, beta, gamma, garbage_col, betaGarbage):
        '''
        construct HMM for motif model
        @params
        negLLMatrix: the negative log likelihood matrix per state
        motif: a list of cluster indices that are the motif. i.e [0,1,0] is A^* B^* A^*
        garbageCol: the garbage column's likelihoods
        betaGarbage: the timesteps where you should add beta to the transition
        '''
        self.betaGarbage = betaGarbage
        self.beta = beta
        self.motif = motif
        self.gamma = gamma
        adjacencyMatrix = self.createAdjacencyMatrix(motif)
        negLLMatrix, initDistribution = self.createNegLLMatrix(
            negLLMatrix, motif, garbage_col)
        HMM.__init__(self, adjacencyMatrix, negLLMatrix, initDistribution)

    def GetAdjacencyMatrix(self, ts):
        ''' override super to take into account betaGarbage'''
        if ts not in self.betaGarbage:
            return self.adjacencyMatrix
        else:
            result = self.adjacencyMatrix.copy()
            result[0,0] = self.beta # add a beta penalty
            return result


    def createAdjacencyMatrix(self, motif):
        '''
        Create a circular adjacency matrix that has 1 garbage state
        and then the motif states
        '''
        num_motif_states = len(motif) + 1
        adj = np.full((num_motif_states, num_motif_states), np.infty)
        np.fill_diagonal(adj, 0)  # allow transition to the same state
        np.fill_diagonal(adj[:, 1:], self.beta)  # allow transition to the next state
        adj[-1, 0] = self.beta  # allow last state to go back to garbage
        adj[-1, 1] = self.beta # allow last state to go back to start state
        return adj

    def createNegLLMatrix(self, negLLMatrix, motif, garbage_col):
        '''
        Create the actual nxeg ll matrix for the motif by extracting out the 
        relative likelihood cols and adding a garbage cols (which is col 0)
        '''
        n, _ = np.shape(negLLMatrix)

        # compute the likelihoods assigned to the garbage value. TODO: figure this out
        # (1) discounted best value
        # (2) average out best values and discount, commented out below
        # (3) discounted original value
        # bestVal = self.gamma*np.mean(np.exp(-1*np.min(negLLMatrix, axis=1)))
        # garbageValue = -1*np.log(bestVal)
        # garbageCol = np.array([[garbageValue]*n]).T
        
        # bestVals = np.min(negLLMatrix, axis=1) - np.log(self.gamma)
        # garbageCol = np.reshape(bestVals, (n, 1))
        # origVals = negLLMatrix[range(n), original_assign] - np.log(self.gamma)
        # origVals = negLLMatrix[np.ix_(np.arange(n), original_assign)] - np.log(self.gamma)
        # garbageCol = np.reshape(origVals, (n, 1))

        garbageCol = np.reshape(garbage_col, (n,1))

        negLLMatrix = np.take(negLLMatrix, motif, axis=1)
        # initial distribution should allow either the garbage or the first state
        initDistribution = np.full(len(motif) + 1, np.infty)
        initDistribution[0:2] = 0
        finalResult = np.c_[garbageCol, negLLMatrix]
        return finalResult, initDistribution

    def SolveAndReturn(self):
        '''
        Solve the HMM and return the most likely sequence

        Returns
        --------
        value1: assignments of each point to a motif state. -1 for garbage
        value2: for each motif instance, tuple of the form (splitpoints, neg log likelihood score)
        '''
        self.solveHMM()
        # get score for garbage and for last motif
        lastStage = len(self.motif)
        garbageScore = self.getEndingScore(0)
        endMotifScore = self.getEndingScore(lastStage)
        if garbageScore < endMotifScore:
            # we ended on garbage
            result = self.GenerateSequenceFromBackPointer(ending=0)
        else:
            result = self.GenerateSequenceFromBackPointer(ending=lastStage)
        result = np.array(result) - 1  # offset so that garbage is -1
        return result, self.SplitAndGetLikelihood(result)

    def SplitAndGetLikelihood(self, sequence):
        # motif list: (index, change, ..., end): likelihood
        motif_length = len(self.motif)
        motifs = []
        currSeq = -1
        currMotif = []
        # get motif changepoints
        for i, c in enumerate(sequence):
            if c != currSeq:
                # if this is the end of the motif, append end point
                if c < currSeq:
                    currMotif.append(i-1)
                    motifs.append(tuple(currMotif))
                    currMotif = []
                    if c != -1:  # straight beginning of new motif
                        currMotif.append(i)
                else:
                    currMotif.append(i)
                currSeq = c
        if c != -1:
            currMotif.append(i)
            motifs.append(tuple(currMotif))

        results = []
        for m in motifs:
            start = m[0]
            end = m[-1]
            ll_end = self.viterbiGrid[end, motif_length-1]
            ll_prev = 0
            if start != 0:
                prev = self.backPointers[start, 1]  # get the start
                ll_prev = self.viterbiGrid[start-1, prev]
            likelihood = ll_end - ll_prev
            results.append((m, likelihood))
        return results
