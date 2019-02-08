#expecting a results.md
#would make the graph for speedup and efficiency vs number of threads
import matplotlib.pyplot as plt
f = open("results_nthreads.md")
seq = 0
pthread = []
omp = []
ctr=0
MODE = 0 #0 for seq, 1 for pthread, 2 for omp
for line in f.readlines():
	if line == "":
		continue
	elif line.strip()=="#seq":
		MODE = 0
	elif line.strip()=="#pthread":
		MODE = 1
	elif line.strip()=="#omp":
		MODE = 2

	line = line.split()
	if line[0][0]=='#':
		continue
	if MODE==0:
		seq = float(line[0])
	elif MODE==1:
		pthread.append(seq/float(line[0]))
	else:
		omp.append(seq/float(line[0]))
	ctr+=1

#pthreads
plt.plot(range(1,len(pthread)+1), pthread, marker=".")
plt.grid()
plt.title("Pthreads: Speedup vs Number of threads")
plt.xlabel("Number of threads")
plt.ylabel("Speedup")
plt.show()

for i in range(len(pthread)):
	pthread[i] /= i+1
plt.plot(range(1,len(pthread)+1), pthread, marker=".")
plt.grid()
plt.title("Pthreads: Efficiency vs Number of threads")
plt.xlabel("Number of threads")
plt.ylabel("Efficiency")
plt.show()

#omp
plt.plot(range(1,len(omp)+1), omp, marker=".")
plt.grid()
plt.title("OMP: Speedup vs Number of threads")
plt.xlabel("Number of threads")
plt.ylabel("Speedup")
plt.show()

for i in range(len(omp)):
	omp[i] /= i+1
plt.plot(range(1,len(omp)+1), omp, marker=".")
plt.grid()
plt.title("OMP: Efficiency vs Number of threads")
plt.xlabel("Number of threads")
plt.ylabel("Efficiency")
plt.show()