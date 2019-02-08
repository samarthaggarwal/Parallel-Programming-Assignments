#expecting a results_problemsize.md
#would make the graph for speedup and efficiency vs number of threads
import matplotlib.pyplot as plt
f = open("results_problemsize.md")
seq = []
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
		ctr=0
		MODE = 1
	elif line.strip()=="#omp":
		ctr=0
		MODE = 2

	line = line.split()
	if line[0][0]=='#':
		continue
	if MODE==0:
		seq.append(float(line[0]))
		# seq = float(line[0])
	elif MODE==1:
		pthread.append(seq[ctr]/float(line[0]))
	else:
		omp.append(seq[ctr]/float(line[0]))
	ctr+=1

#pthreads
plt.plot(range(1000,1000*len(pthread)+1,1000), pthread, marker=".")
plt.grid()
plt.title("Pthreads: Speedup vs Problem Size")
plt.xlabel("Problem size")
plt.ylabel("Speedup")
plt.show()

for i in range(len(pthread)):
	pthread[i] /= 4
plt.plot(range(1000,1000*len(pthread)+1,1000), pthread, marker=".")
plt.grid()
plt.title("Pthreads: Efficiency vs Problem Size")
plt.xlabel("Problem size")
plt.ylabel("Efficiency")
plt.show()

#omp
plt.plot(range(1000,1000*len(omp)+1,1000), omp, marker=".")
plt.grid()
plt.title("OMP: Speedup vs Problem Size")
plt.xlabel("Problem size(N)")
plt.ylabel("Speedup")
plt.show()

for i in range(len(omp)):
	omp[i] /= 4
plt.plot(range(1000,1000*len(omp)+1,1000), omp, marker=".")
plt.grid()
plt.title("OMP: Efficiency vs Problem Size")
plt.xlabel("Problem size(N)")
plt.ylabel("Efficiency")
plt.show()