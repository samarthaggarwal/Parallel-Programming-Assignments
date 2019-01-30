all:
	g++ kmeans2.cpp -o kmeans.out -fopenmp
	./kmeans.out < in.txt

pthread:
	g++ kmeans_pthread.cpp -o kmeans_pthread.out -pthread -fopenmp
	./kmeans_pthread.out < in.txt

openmp:
	g++ kmeans_openmp.cpp -o kmeans_openmp.out -fopenmp
	./kmeans_openmp.out < in3.txt

clean:
	rm *.out

randomInput:
	g++ randomInput.cpp
	./a.out > in1.txt

