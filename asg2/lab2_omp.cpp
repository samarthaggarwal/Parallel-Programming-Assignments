// #include <malloc.h>
#include<stdlib.h>
#include <omp.h>
#include<stdio.h>
#include<iostream>
#include<algorithm>
#include<cmath>
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

void print_vector(float *a, int dim){
	for(int i=0;i<dim;i++)
		printf("%f\t", a[i]);
	printf("\n");
}

float dot_product(float *a, float *b, int dim){
	// cout<<"called dot_product bw\n";
	// print_vector(a,dim);
	// print_vector(b,dim);
	float sum=0;
	for(int i=0;i<dim;i++){
		sum+=a[i]*b[i];
	}

	// printf("dot=%f\n", sum);
	return sum;
}

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
	float a[N][N];

// making d and dTranspose
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			d[i][j] = D[M*i+j];
			dTranspose[j][i] = D[M*i+j];
		}
	}
/*
// print d
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			// printf("%d %d \t\t", i,j);
			printf("%.0f\t", d[i][j] );
		}
		printf("\n");
	}
	printf("\n");

// print dTranspose
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			// printf("%d %d \t\t", i,j);
			printf("%.0f\t", dTranspose[i][j] );
		}
		printf("\n");
	}
	printf("\n");

// calc dTd
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			a[i][j]=0;
			for(int k=0;k<M;k++){
				a[i][j] += dTranspose[i][k]*d[k][j];
			}
		}
	}

// print dTd
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			// printf("%d %d \t\t", i,j);
			printf("%.0f\t", a[i][j] );
		}
		printf("\n");
	}
	printf("\n");
*/
// QR decomposition
	cout<<"A\n";
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			// printf("%d %d \t\t", i,j);
			printf("%.0f\t", d[i][j] );
		}
		printf("\n");
	}

	cout<<"transposing A\n";
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			a[i][j]=dTranspose[i][j];
		}
	}

	cout<<"A transpose\n";
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			// printf("%d %d \t\t", i,j);
			printf("%.0f\t", a[i][j] );
		}
		printf("\n");
	}
	printf("\n");
	float q[N][N];
	float q2[N][N];
	float r[N][N];
	float temp[N];
	float dot,normSquared,norm;

	for(int i=0;i<N;i++){
		// cout<<"i="<<i<<endl;
		for(int j=0;j<N;j++){
			q[i][j]=a[i][j];
		}

		// cout<<"init q[i] as ";
		// print_vector(q[i],N);
		
		for(int k=0;k<i;k++){
			dot=dot_product(a[i],q[k],N);
			normSquared=dot_product(q[k],q[k],N);
			// cout<<"k="<<k<<"\tdot="<<dot<<"\tnormSq="<<normSquared<<endl;
			for(int j=0;j<N;j++){
				// cout<<"step "<<q[i][j]<<"\t"<<(a[i][j]*dot)/normSquared << "\t"<<q[i][j] - (a[i][j]*dot)/normSquared<<endl;
				q[i][j] -= (q[k][j]*dot)/normSquared;
			}
		}

		// cout<<"updated q[i] as ";
		// print_vector(q[i],N);

		// norm=sqrt(dot_product(q[i],q[i],N));
		// for(int j=0;j<N;j++){
		// 	q[i][j]/=norm;
		// }
	}
	
	cout<<"U transpose\n";
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			// printf("%d %d \t\t", i,j);
			printf("%f\t", q[i][j] );
		}
		printf("\n");
	}
	printf("\n");

	cerr<<"normalising\n";
	for(int i=0;i<N;i++){
		norm=sqrt(dot_product(q[i],q[i],N));
		for(int j=0;j<N;j++){
			q[i][j]/=norm;
		}
	}

	cerr<<"transposing\n";
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			q2[i][j]=q[j][i];
		}
	}

	cout<<"Q\n";
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			// printf("%d %d \t\t", i,j);
			printf("%f\t", q2[i][j] );
		}
		printf("\n");
	}
	printf("\n");	

	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			if(j<i){
				r[i][j]=0;
				continue;
			} else{
				r[i][j]=dot_product(q[i],a[j],N);
			}
		}
	}

	cout<<"R\n";
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			// printf("%d %d \t\t", i,j);
			printf("%f\t", r[i][j] );
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
