all:
	g++ -fopenmp src/main_sequential.c src/lab1_sequential.cpp src/lab1_io.c
	if [[ $1 -eq 5 ]]; then
		./a.out 3 testing/test5.txt outputs/points.txt outputs/centres.txt
	elif [[ $1 -eq 50 ]]; then
		./a.out 4 testing/test50.txt outputs/points.txt outputs/centres.txt
	elif [[ $1 -eq 1l ]]; then
		./a.out 3 testing/test1l.txt outputs/points.txt outputs/centres.txt
	fi
	# ./a.out 3 testing/test5.txt outputs/points.txt outputs/centres.txt
	python3 visualise.py outputs/points.txt

pthread:
	g++ -fopenmp -pthread src/main_pthread.c src/lab1_pthread.cpp src/lab1_io.c
	./a.out 3 4 testing/test5.txt outputs/points.txt outputs/centres.txt

openmp:
	g++ -fopenmp -pthread src/main_omp.c src/lab1_omp.cpp src/lab1_io.c
	./a.out 4 4 testing/test50.txt outputs/points.txt outputs/centres.txt
