# testing on 50K test case (4 clusters)

# numThreads=4
# name="dataset_"
# suffix="_10.txt"

# for (( i = 10000; i <= 100000; i+=10000 )); do
# 	echo $name$i$suffix
# done

echo "testing seq"
g++ -fopenmp src/main_sequential.c src/lab1_sequential.cpp src/lab1_io.c
./a.out 10 testing/dataset_10000_10.txt outputs/points1.txt outputs/centres1.txt
./a.out 10 testing/dataset_20000_10.txt outputs/points2.txt outputs/centres2.txt
./a.out 10 testing/dataset_30000_10.txt outputs/points3.txt outputs/centres3.txt
./a.out 10 testing/dataset_40000_10.txt outputs/points4.txt outputs/centres4.txt
./a.out 10 testing/dataset_50000_10.txt outputs/points5.txt outputs/centres5.txt
./a.out 10 testing/dataset_60000_10.txt outputs/points6.txt outputs/centres6.txt
./a.out 10 testing/dataset_70000_10.txt outputs/points7.txt outputs/centres7.txt
./a.out 10 testing/dataset_80000_10.txt outputs/points8.txt outputs/centres8.txt
./a.out 10 testing/dataset_90000_10.txt outputs/points9.txt outputs/centres9.txt
./a.out 10 testing/dataset_100000_10.txt outputs/points10.txt outputs/centres10.txt
echo "\n===============\n"

echo "testing pth"
g++ -fopenmp -pthread src/main_pthread.c src/lab1_pthread.cpp src/lab1_io.c
./a.out 10 4 testing/dataset_10000_10.txt outputs/points1.txt outputs/centres1.txt
./a.out 10 4 testing/dataset_20000_10.txt outputs/points2.txt outputs/centres2.txt
./a.out 10 4 testing/dataset_30000_10.txt outputs/points3.txt outputs/centres3.txt
./a.out 10 4 testing/dataset_40000_10.txt outputs/points4.txt outputs/centres4.txt
./a.out 10 4 testing/dataset_50000_10.txt outputs/points5.txt outputs/centres5.txt
./a.out 10 4 testing/dataset_60000_10.txt outputs/points6.txt outputs/centres6.txt
./a.out 10 4 testing/dataset_70000_10.txt outputs/points7.txt outputs/centres7.txt
./a.out 10 4 testing/dataset_80000_10.txt outputs/points8.txt outputs/centres8.txt
./a.out 10 4 testing/dataset_90000_10.txt outputs/points9.txt outputs/centres9.txt
./a.out 10 4 testing/dataset_100000_10.txt outputs/points10.txt outputs/centres10.txt
echo "\n===============\n"

echo "testing omp"
g++ -fopenmp -pthread src/main_omp.c src/lab1_omp.cpp src/lab1_io.c
./a.out 10 4 testing/dataset_10000_10.txt outputs/points1.txt outputs/centres1.txt
./a.out 10 4 testing/dataset_20000_10.txt outputs/points2.txt outputs/centres2.txt
./a.out 10 4 testing/dataset_30000_10.txt outputs/points3.txt outputs/centres3.txt
./a.out 10 4 testing/dataset_40000_10.txt outputs/points4.txt outputs/centres4.txt
./a.out 10 4 testing/dataset_50000_10.txt outputs/points5.txt outputs/centres5.txt
./a.out 10 4 testing/dataset_60000_10.txt outputs/points6.txt outputs/centres6.txt
./a.out 10 4 testing/dataset_70000_10.txt outputs/points7.txt outputs/centres7.txt
./a.out 10 4 testing/dataset_80000_10.txt outputs/points8.txt outputs/centres8.txt
./a.out 10 4 testing/dataset_90000_10.txt outputs/points9.txt outputs/centres9.txt
./a.out 10 4 testing/dataset_100000_10.txt outputs/points10.txt outputs/centres10.txt
echo "\n===============\n"


