# python3 dataset.py
# mv dataset_* testing/testr.txt
# g++ -fopenmp src/main_sequential.c src/lab1_sequential.cpp src/lab1_io.c
# ./a.out 3 testing/testr.txt outputs/points.txt outputs/centres.txt
# python3 visualise.py outputs/points.txt

if [[ $1 = "seq" ]]; then
	g++ -fopenmp src/main_sequential.c src/lab1_sequential.cpp src/lab1_io.c
	if [[ $2 -eq 5 ]]; then
		./a.out 3 testing/test5.txt outputs/points.txt outputs/centres.txt
	elif [[ $2 -eq 50 ]]; then
		./a.out 4 testing/test50.txt outputs/points.txt outputs/centres.txt
	elif [[ $2 = "1l" ]]; then
		./a.out 7 testing/test1l.txt outputs/points.txt outputs/centres.txt
	fi

elif [[ $1 = "pth" ]]; then
	g++ -fopenmp -pthread src/main_pthread.c src/lab1_pthread.cpp src/lab1_io.c
	if [[ $2 -eq 5 ]]; then
		./a.out 3 $3 testing/test5.txt outputs/points.txt outputs/centres.txt
	elif [[ $2 -eq 50 ]]; then
		./a.out 4 $3 testing/test50.txt outputs/points.txt outputs/centres.txt
	elif [[ $2 = "1l" ]]; then
		./a.out 10 $3 testing/test1l.txt outputs/points.txt outputs/centres.txt
	fi

elif [[ $1 = "omp" ]]; then
	if [[ $2 -eq 5 ]]; then
		./a.out 3 testing/test5.txt outputs/points.txt outputs/centres.txt
	elif [[ $2 -eq 50 ]]; then
		./a.out 4 testing/test50.txt outputs/points.txt outputs/centres.txt
	elif [[ $2 = "1l" ]]; then
		./a.out 7 testing/test1l.txt outputs/points.txt outputs/centres.txt
	fi

fi

# ./a.out 3 testing/test5.txt outputs/points.txt outputs/centres.txt
python3 visualise.py outputs/points.txt

