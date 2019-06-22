import numpy as np

size = 1000
num = "00"
results = []

def trial():
	s = np.random.randint(0,2,size).tolist()
	s = [str(i) for i in s]
	s = "".join(s)
	results.append(countNonOverlapping(s, num))

def countNonOverlapping(s, p):
	pos = -1 
	count = 0
	while True:
		pos = s.find(p, pos + 1)
		if pos == -1: break
		count += 1
	return count

for t in range(100000):
	if t % 100 == 0: print(t)
	trial()
print(np.median(results))
