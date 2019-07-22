# mpic++ -lm main_mpi.c lab4_mpi.cpp lab4_io.c -o ppm
# # mpirun -np 1 ./ppm ./testcase/test_pdf
# mpirun -np 4 ./ppm ./testcase/man_test
# # mpirun -np 4 ./ppm ./testcase/testcase_10000_10_aval
# # mpirun -np 1 ./ppm ./testcase/testcase_1000000_100

# mpirun -np 1 ./ppm ./testcase/test_pdf
# mpirun -np 4 ./ppm ./testcase/testcase_10000_10_aval
# mpirun --oversubscribe -np 7 ./ppm ./testcase/testcase_1000000_100

mpic++ -lm main_mpi.c lab4_mpi.cpp lab4_io.c -o sam
mpic++ -lm main_mpi.c lab4_mpi_ati.c lab4_io.c -o ati

mpirun -np $2 ./sam $1 > sam.txt
mpirun -np $2 ./ati $1 > ati.txt

echo "sam vs ati"
diff sam.txt ati.txt
