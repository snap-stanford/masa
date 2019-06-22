import matplotlib.pyplot as plt 



y = [107.84717392921448, 189.18333387374878, 426.1840078830719, 909.5116035938263, 1501.5912849903107, 2479.0617039203644, 3251.2408401966095, 3696.750289916992]
x = [1000, 1500, 2000, 3000, 4000, 5000, 6000, 7000]
x = [i*15*10 for i in x]
plt.figure(1, figsize=(7, 4))
plt.plot(x,y, "-bs", markersize=4)
plt.xlabel("T")
plt.ylabel("Per-iteration runtime (s)")
plt.rcParams.update({'font.size': 14})
plt.savefig('scalability.eps', format='eps', bbox_inches='tight', dpi=1000)
plt.show()