# K-means-clustering
COL380 - lab1

## Sub :
1. add readme containing flags reqd for compilation


### Opti Seq:
- [ ] pthreads - done variable (array) to sync instead of spawning and killing threads


### Report:
1. Why +1 in lengthPerThread - load balancing
2. Why not parallelised recompute_means -> need locks, will be worse than sequential
	- lock inside for loop - worse than sequential code
3. Testing
	- by plotting graphs (sankalan's script) and checking visualisations
	- by finding serial fraction of code from diff speedups. it should match
	


