all:
	g++ kmeans2.cpp -o kmeans.out
	./kmeans.out < in.txt

pthread:
	g++ kmeans_pthread.cpp -o kmeans_pthread.out -pthread -fopenmp
	./kmeans_pthread.out < in.txt

openmp:
	g++ kmeans_openmp.cpp -o kmeans_openmp
	./kmeans_openmp.out < in.txt

clean:
	rm *.out

randomInput:
	g++ randomInput.cpp
	./a.out > in1.txt

