import numpy as np
import matplotlib
# matplotlib.use('agg')
import matplotlib.pyplot as plt
from analyze_synthetic import getValidMappings 
print("Started")
GAMMAS = ["0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "0.99"]

def getScores(directory):
    correctname = "%s/correct.out" % directory
    scores = [[],[],[]]
    for g in GAMMAS:
        print(g)
        print(correctname)
        s = getValidMappings(correctname,  "%s/0.2/%s/assign.out" % (directory, g))
        for i in range(len(s)): scores[i].append(s[i])
    ticcscores = getValidMappings(correctname,  "%s/0.2/old/assign.out" % (directory))
    return scores, ticcscores

def saveScores(directory, output):
    all_scores = getScores(directory)
    np.savetxt("%s/%s" % (output, "gammaRobust.out"), np.array(all_scores), delimiter=',')

def loadAndPlotScores(directory):
    scoreValues, ticcscores = getScores(directory)
    print(scoreValues, ticcscores)
    labels = [str(v) for v in GAMMAS]
    titles = ["$\gamma$ Robustness: Weighted Macro F1 Score on Motif Clusters", 
        "Weighted Macro F1 Score on All Clusters",
        "Accuracy Score on Motif Segments"]
    for i in range(0,1):
        plt.figure(i, figsize=(7, 3))
        plt.plot(labels, scoreValues[i], "-bs", label="CASC")
        ticcLine = [ticcscores[i] for _ in range(len(labels))]
        plt.plot(labels, ticcLine, "--", label="TICC", color="orange")
        plt.xlabel("$\gamma$")
        plt.ylim(ymin=0.5, ymax=1)
        # plt.title(titles[i])
        if i == 2:
            plt.ylabel("Accuracy score on motif segments vs Noise")
        else:
            plt.ylabel("Weighted F1 Score")
        # plt.legend()
        plt.legend(loc='lower center', ncol=2, fancybox=False, edgecolor="black")
        plt.rcParams.update({'font.size': 14})
        plt.savefig('robustness.eps', format='eps', bbox_inches='tight', dpi=1000)

    plt.show()


# saveScores("ordered_synthetic_retry2", "final_scores_retry2")
loadAndPlotScores("ordered_synthetic_allperturbs/")


