echo "running sequential\n"
g++ -fopenmp src/main_sequential.c src/lab1_sequential.cpp src/lab1_io.c
./a.out 10 testing/test5.txt outputs/points.txt outputs/centres.txt
echo "========"
./a.out 10 testing/test50.txt outputs/points.txt outputs/centres.txt
echo "========"
./a.out 10 testing/test10l.txt outputs/points.txt outputs/centres.txt
echo "========"
