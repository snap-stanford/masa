import matplotlib.pyplot as plt

"""Accuracy plots"""

labels = [0.4,0.5,0.6,0.7,0.9,1]

#F1 motif
plt.figure(1)
ax = plt.subplot(111)
TICCmotif = [0.857,0.786,0.783,0.812,0.804,0.839]
CASCmotif = [0.894,0.854,0.842,0.851,0.829,0.856]
TICCTotal = [0.802,0.673,0.671,0.691,0.684,0.712]
CASCTotal = [0.825,0.713,0.708,0.715,0.702,0.728]
plt.plot(labels, CASCmotif, marker='o', color='r', label="CASC, scored on motif clusters")
plt.plot(labels, TICCmotif, marker='*', color='b', label="TICC, scored on motif clusters")
plt.plot(labels, CASCTotal, marker='o', color='r', linestyle='dashed', label="CASC, scored on all clusters")
plt.plot(labels, TICCTotal, marker='*', color='b', linestyle='dashed', label="TICC, scored on all clusters")
plt.title("CASC vs TICC F1 score")
# plt.legend()

box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.15,
                 box.width, box.height * 0.85])

ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.14),
          fancybox=True, shadow=True, ncol=2)
plt.xlabel("Fraction of noisy segments")
plt.ylabel("Weighted F1 score")

# Accuracy
plt.figure(2)
ax = plt.subplot(111)
TICCAccuracy = [0.778, 0.661, 0.654, 0.691, 0.659, 0.700]
CASCAccuracy = [0.901, 0.892, 0.875, 0.856, 0.830, 0.793]
plt.plot(labels, CASCAccuracy, marker='o', color='r', label="CASC")
plt.plot(labels, TICCAccuracy, marker='*', color='b', label="TICC")
plt.title("CASC vs TICC Accuracy on motif segments")
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.15,
                 box.width, box.height * 0.85])
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.14),
          fancybox=True, shadow=True, ncol=2)
plt.xlabel("Fraction of noisy segments")
plt.ylabel("Accuracy on motif segments")

plt.show()