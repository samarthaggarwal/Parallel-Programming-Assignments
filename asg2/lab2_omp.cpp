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

	// float d[M][N];
	float **d=(float**)malloc(sizeof(float*)*M);
	for(int i=0;i<M;i++) d[i]=(float*)malloc(sizeof(float)*N);

	// float dTranspose[N][M];
	float **dTranspose=(float**)malloc(sizeof(float*)*N);
	for(int i=0;i<N;i++) dTranspose[i]=(float*)malloc(sizeof(float)*M);

	// float a[N][N];
	float **a=(float**)malloc(sizeof(float*)*N);
	for(int i=0;i<N;i++) a[i]=(float*)malloc(sizeof(float)*N);

// making d and dTranspose
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			d[i][j] = D[N*i+j];
			dTranspose[j][i] = D[N*i+j];
		}
	}

// print d
	// for(int i=0;i<M;i++){
	// 	for(int j=0;j<N;j++){
	// 		// printf("%d %d \t\t", i,j);
	// 		printf("%.0f\t", d[i][j] );
	// 	}
	// 	printf("\n");
	// }
	// printf("\n");

// print dTranspose
	// for(int i=0;i<M;i++){
	// 	for(int j=0;j<N;j++){
	// 		// printf("%d %d \t\t", i,j);
	// 		printf("%.0f\t", dTranspose[i][j] );
	// 	}
	// 	printf("\n");
	// }
	// printf("\n");

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
	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<N;j++){
	// 		// printf("%d %d \t\t", i,j);
	// 		printf("%.0f\t", a[i][j] );
	// 	}
	// 	printf("\n");
	// }
	// printf("\n");

// QR decomposition
	// float q[N][N];
	float **q=(float**)malloc(sizeof(float*)*N);
	for(int i=0;i<N;i++) q[i]=(float*)malloc(sizeof(float)*N);

	// float q2[N][N];
	float **q2=(float**)malloc(sizeof(float*)*N);
	for(int i=0;i<N;i++) q2[i]=(float*)malloc(sizeof(float)*N);
	
	// float r[N][N];
	float **r=(float**)malloc(sizeof(float*)*N);
	for(int i=0;i<N;i++) r[i]=(float*)malloc(sizeof(float)*N);

	// float e[N][N];
	float **e=(float**)malloc(sizeof(float*)*N);
	for(int i=0;i<N;i++) e[i]=(float*)malloc(sizeof(float)*N);

	// float eUpdated[N][N];
	float **eUpdated=(float**)malloc(sizeof(float*)*N);
	for(int i=0;i<N;i++) eUpdated[i]=(float*)malloc(sizeof(float)*N);

	// float temp[N];
	float *temp=(float*)malloc(sizeof(float*)*N);

	// float sigma[N][N];
	float **sigma=(float**)malloc(sizeof(float*)*N);
	for(int i=0;i<N;i++) sigma[i]=(float*)malloc(sizeof(float)*N);

	// float sigmaInv[N][N];
	float **sigmaInv=(float**)malloc(sizeof(float*)*N);
	for(int i=0;i<N;i++) sigmaInv[i]=(float*)malloc(sizeof(float)*N);

	// float u[M][N];
	float **u=(float**)malloc(sizeof(float*)*M);
	for(int i=0;i<M;i++) u[i]=(float*)malloc(sizeof(float)*N);

	float dot,normSquared,norm, maxDiff;

	// setting a = d
	// cout<<"A\n";
	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<N;j++){
	// 		// printf("%d %d \t\t", i,j);
	// 		printf("%.0f\t", d[i][j] );
	// 	}
	// 	printf("\n");
	// }
	// printf("\n");

	// cout<<"transposing A\n";
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			// a[i][j]=dTranspose[i][j];
			e[i][j]=0;
		}
		e[i][i]=1;
	}

	// cout<<"A transpose\n";
	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<N;j++){
	// 		// printf("%d %d \t\t", i,j);
	// 		printf("%.0f\t", a[i][j] );
	// 	}
	// 	printf("\n");
	// }

	int numIter=0;
	do{

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

			norm=sqrt(dot_product(q[i],q[i],N));
			for(int j=0;j<N;j++){
				q[i][j]/=norm;
			}
		}
		
		// cout<<"U transpose\n";
		// for(int i=0;i<N;i++){
		// 	for(int j=0;j<N;j++){
		// 		// printf("%d %d \t\t", i,j);
		// 		printf("%f\t", q[i][j] );
		// 	}
		// 	printf("\n");
		// }
		// printf("\n");

		// cerr<<"normalising\n";
		// for(int i=0;i<N;i++){
		// 	norm=sqrt(dot_product(q[i],q[i],N));
		// 	for(int j=0;j<N;j++){
		// 		q[i][j]/=norm;
		// 	}
		// }


		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				if(j>i){
					r[i][j]=0;
					continue;
				} else{
					r[i][j]=dot_product(q[j],a[i],N);
				}
			}
		}

		cout<<"R transpose\n";
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				// printf("%d %d \t\t", i,j);
				printf("%f\t", r[i][j] );
			}
			printf("\n");
		}
		printf("\n");

		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				a[i][j]=0;
				eUpdated[i][j]=0;
				for(int k=0;k<N;k++){
					a[i][j]+=q[i][k]*r[k][j];
					eUpdated[i][j]+=q[i][k]*e[k][j];
				}
			}
		}

		maxDiff=0;
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				if(fabs(e[i][j]-eUpdated[i][j]) > maxDiff)
					maxDiff=fabs(e[i][j]-eUpdated[i][j]);
				e[i][j]=eUpdated[i][j];
			}
		}

		numIter++;
		cout<<"iter = "<<numIter<<endl;
	}while(numIter < 2);

	cout<<"numIter = "<<numIter<<endl;

	// cerr<<"transposing\n";
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

// eUpdated has eTranspose now
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			eUpdated[i][j]=e[j][i];
		}
	}

	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			sigma[i][j]=0;
			sigmaInv[i][j]=0;
		}
		sigma[i][i]=sqrt(a[i][i]);
		sigmaInv[i][i]=1/sigma[i][i];
	}

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			u[i][j]=0;
			for(int k=0;k<N;k++){
				u[i][j]+=d[i][k]*eUpdated[k][j];
			}
			u[i][j]*=sigmaInv[j][j];
		}
	}

	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			(*U)[N*i+j]=eUpdated[i][j];
		}
	}

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			(*V_T)[N*i+j]=u[i][j];
		}
	}
	for(int i=M;i<N;i++){
		for(int j=0;j<N;j++){
			(*V_T)[N*i+j]=0;
		}
	}
	
	for(int i=0;i<N;i++){
		*(*SIGMA +i)=sigma[i][i];
	}

}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    printf("called pca\n");

    float* sigmaSorted = (float*)malloc(sizeof(float)*N);
    float sum=0, ret=retention/100.0;

    for(int i=0;i<N;i++){
    	sigmaSorted[i]=SIGMA[i];
    	sum+=sigmaSorted[i];
    }
    for(int i=0;i<N;i++){
    	sigmaSorted[i]/=sum;
    }
    sort(sigmaSorted,sigmaSorted + N);

    for(int i=0;i<N;i++){
    	ret-=sigmaSorted[i];
    	if(ret<=0){
    		*K=i+1;
    		break;
    	}
    }

    // float w[N][K];
	float **w=(float**)malloc(sizeof(float*)*N);
	for(int i=0;i<N;i++) w[i]=(float*)malloc(sizeof(float)* *K);

	int j;
    for(int i=0;i<*K;i++){
    	for(j=0;j<N;j++){
    		if(fabs(sigmaSorted[i]-SIGMA[j])<0.0001 )
    			break;
    	}
    	for(int k=0;k<N;k++){
    		w[i][k]=U[N*j+k];
    	}
    }

    *D_HAT = (float*)malloc(sizeof(float)*M* *K);
    for(int i=0;i<M;i++){
    	for(int j=0;j<*K;j++){
    		(*D_HAT)[*K *i + j] =0;
    		for(int l=0;l<N;l++){
	    		(*D_HAT)[*K *i+j] += D[N*i + l]*w[l][j];
    		}
    	}
    }

    // for(int i=0;i<M;i++){
    // 	for(int j=0;j<*K;j++){
    // 		printf("%.6f\t", (*D_HAT)[*K *i+j]);
    // 	}
    // 	printf("\n");
    // }
}
