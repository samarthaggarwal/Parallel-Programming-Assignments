# testing on 50K test case (4 clusters)

numThreads=1

echo "testing seq"
g++ -fopenmp src/main_sequential.c src/lab1_sequential.cpp src/lab1_io.c
./a.out 4 testing/test50.txt outputs/points.txt outputs/centres.txt

echo "testing pth"
g++ -fopenmp -pthread src/main_pthread.c src/lab1_pthread.cpp src/lab1_io.c
for (( i = 1; i <= 10; i++ )); do
	./a.out 4 $i testing/test50.txt outputs/points.txt outputs/centres.txt
done

echo "testing omp"
g++ -fopenmp -pthread src/main_omp.c src/lab1_omp.cpp src/lab1_io.c
for (( i = 1; i <= 10; i++ )); do
	./a.out 4 $i testing/test50.txt outputs/points.txt outputs/centres.txt
done

