all:
	g++ -fopenmp src/main_sequential.c src/lab1_sequential.cpp src/lab1_io.c
	./a.out 10 testing/test4.txt outputs/out1.txt outputs/out2.txt

pthread:
	g++ -fopenmp -pthread src/main_pthread.c src/lab1_pthread.cpp src/lab1_io.c
	./a.out 10 4 testing/test3.txt outputs/out1.txt outputs/out2.txt

openmp:
	g++ -fopenmp -pthread src/main_omp.c src/lab1_omp.cpp src/lab1_io.c
	./a.out 10 1 testing/test3.txt outputs/out1.txt outputs/out2.txt
