// #include <malloc.h>
#include<stdlib.h>
#include <omp.h>
#include<stdio.h>
#include<iostream>
#include<algorithm>
using namespace std;

// a is p x q , b is q x r, res is p x r
// void matmul(float **a,float **b, float **res, int p, int q, int r){
// 	for (int i=0;i<p;i++){
// 		for (int j=0;j<r;j++){
// 			res[i][j]=0;
// 			for(int k=0;k<q;k++){
// 				res[i][j] += a[i][k]*b[k][j];
// 			}
// 		}
// 	}
// }

// void print(void *mat, int a, int b){
// 	printf("called print\n");
// 	float *m[b] = (float *[m])mat;
// 	for(int i=0;i<a;i++){
// 		for(int j=0;j<b;j++){
// 			printf("%d %d \t\t", i,j);
// 			// printf("%f\n", *(*(m+i)+j) );
// 		}
// 		printf("\n");
// 	}
// }

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{
	printf("called svd\n");

	float d[M][N];
	float dTranspose[N][M];
	float dTd[N][N];

// making d and dTranspose
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			d[i][j] = D[M*i+j];
			dTranspose[j][i] = D[M*i+j];
		}
	}

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			// printf("%d %d \t\t", i,j);
			printf("%.0f\t", d[i][j] );
		}
		printf("\n");
	}
	printf("\n");

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			// printf("%d %d \t\t", i,j);
			printf("%.0f\t", dTranspose[i][j] );
		}
		printf("\n");
	}
	printf("\n");

	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			dTd[i][j]=0;
			for(int k=0;k<M;k++){
				dTd[i][j] += dTranspose[i][k]*d[k][j];
			}
		}
	}

	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			// printf("%d %d \t\t", i,j);
			printf("%.0f\t", dTd[i][j] );
		}
		printf("\n");
	}
	printf("\n");
	
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    printf("called pca\n");
}
