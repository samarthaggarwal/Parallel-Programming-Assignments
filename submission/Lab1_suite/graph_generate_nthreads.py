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

# for p in range(1,10):
# 	print(pthread[p])
# 	print(omp[p])
# 	s = pthread[p-1]
# 	print( "f = " + str(( (1/s)-(1/(p+1)) )/( 1 - (1/(p+1)) ) ) )
# 	s = omp[p-1]
# 	print( "f = " + str(( (1/s)-(1/(p+1)) )/( 1 - (1/(p+1)) ) ) )


#pthreads
plt.plot(range(1,len(pthread)+1), pthread, marker=".")
plt.grid()
plt.title("Pthreads: Speedup vs Number of threads")
plt.xlabel("Number of threads")
plt.ylabel("Speedup")
# plt.show()
plt.savefig('outputs/speedup_vs_threads_pthread.png')
plt.close()


for i in range(len(pthread)):
	pthread[i] /= i+1
plt.plot(range(1,len(pthread)+1), pthread, marker=".")
plt.grid()
plt.title("Pthreads: Efficiency vs Number of threads")
plt.xlabel("Number of threads")
plt.ylabel("Efficiency")
# plt.show()
plt.savefig('outputs/efficiency_vs_threads_pthread.png')
plt.close()

#omp
plt.plot(range(1,len(omp)+1), omp, marker=".")
plt.grid()
plt.title("OMP: Speedup vs Number of threads")
plt.xlabel("Number of threads")
plt.ylabel("Speedup")
# plt.show()
plt.savefig('outputs/speedup_vs_threads_omp.png')
plt.close()

for i in range(len(omp)):
	omp[i] /= i+1
plt.plot(range(1,len(omp)+1), omp, marker=".")
plt.grid()
plt.title("OMP: Efficiency vs Number of threads")
plt.xlabel("Number of threads")
plt.ylabel("Efficiency")
# plt.show()
plt.savefig('outputs/efficiency_vs_threads_omp.png')
plt.close()

