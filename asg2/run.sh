#!/bin/bash
g++ -fopenmp -lm lab2_io.c lab2_omp.cpp main_omp.c -o pca
# time ./pca testcase/testcase_200_30 80
