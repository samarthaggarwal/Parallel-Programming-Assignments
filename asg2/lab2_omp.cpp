// #include <malloc.h>
#include<stdlib.h>
#include<omp.h>
#include<stdio.h>
#include<iostream>
#include<algorithm>
#include<cmath>
using namespace std;

#define NUM_THREADS 4

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
// 	//printf("called print\n");
// 	float *m[b] = (float *[m])mat;
// 	for(int i=0;i<a;i++){
// 		for(int j=0;j<b;j++){
// 			//printf("%d %d \t\t", i,j);
// 			// //printf("%f\n", *(*(m+i)+j) );
// 		}
// 		//printf("\n");
// 	}
// }

// void print_vector(double *a, int dim){
	// for(int i=0;i<dim;i++)
		// //printf("%f\t", a[i]);
	// //printf("\n");
// }

double dot_product(double *a, double *b, int dim){
	// cout<<"called dot_product bw\n";
	// print_vector(a,dim);
	// print_vector(b,dim);
	double sum=0;
	for(int i=0;i<dim;i++){
		sum+=a[i]*b[i];
	}

	// //printf("dot=%f\n", sum);
	return sum;
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{

	double startTime = omp_get_wtime();

	printf("num_threads = %d\n", NUM_THREADS);
	omp_set_dynamic(0);
	omp_set_num_threads(NUM_THREADS);
	//printf("called svd\n");

	// double d[M][N];
	double **d=(double**)malloc(sizeof(double*)*M);
	for(int i=0;i<M;i++) d[i]=(double*)malloc(sizeof(double)*N);

	// double dTranspose[N][M];
	double **dTranspose=(double**)malloc(sizeof(double*)*N);
	for(int i=0;i<N;i++) dTranspose[i]=(double*)malloc(sizeof(double)*M);

	// double a[N][N];
	double **a=(double**)malloc(sizeof(double*)*N);
	for(int i=0;i<N;i++) a[i]=(double*)malloc(sizeof(double)*N);

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
	// 		// //printf("%d %d \t\t", i,j);
	// 		//printf("%.0f\t", d[i][j] );
	// 	}
	// 	//printf("\n");
	// }
	// //printf("\n");

// print dTranspose
	// for(int i=0;i<M;i++){
	// 	for(int j=0;j<N;j++){
	// 		// //printf("%d %d \t\t", i,j);
	// 		//printf("%.0f\t", dTranspose[i][j] );
	// 	}
	// 	//printf("\n");
	// }
	// //printf("\n");

// calc dTd
	#pragma omp parallel num_threads(NUM_THREADS)
	{
		#pragma omp for
		for (int i=0;i<N;i++){
			for (int j=0;j<N;j++){
				a[i][j]=0;
				for(int k=0;k<M;k++){
					a[i][j] += dTranspose[i][k]*d[k][j];
				}
			}
		}
	}
// print dTd
	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<N;j++){
	// 		// //printf("%d %d \t\t", i,j);
	// 		//printf("%.0f\t", a[i][j] );
	// 	}
	// 	//printf("\n");
	// }
	// //printf("\n");

// QR decomposition
	// double q[N][N];
	double **q=(double**)malloc(sizeof(double*)*N);
	for(int i=0;i<N;i++) q[i]=(double*)malloc(sizeof(double)*N);

	// double q2[N][N];
	// double **q2=(double**)malloc(sizeof(double*)*N);
	// for(int i=0;i<N;i++) q2[i]=(double*)malloc(sizeof(double)*N);
	
	// double r[N][N];
	double **r=(double**)malloc(sizeof(double*)*N);
	for(int i=0;i<N;i++) r[i]=(double*)malloc(sizeof(double)*N);

	// double e[N][N];
	double **e=(double**)malloc(sizeof(double*)*N);
	for(int i=0;i<N;i++) e[i]=(double*)malloc(sizeof(double)*N);

	// double eUpdated[N][N];
	double **eUpdated=(double**)malloc(sizeof(double*)*N);
	for(int i=0;i<N;i++) eUpdated[i]=(double*)malloc(sizeof(double)*N);

	double **aUpdated=(double**)malloc(sizeof(double*)*N);
	for(int i=0;i<N;i++) aUpdated[i]=(double*)malloc(sizeof(double)*N);

	// double temp[N];
	double *temp=(double*)malloc(sizeof(double*)*N);

	// double aPrev[N];
	double *aPrev=(double*)malloc(sizeof(double*)*N);

	// double sigma[N][N];
	double **sigma=(double**)malloc(sizeof(double*)*N);
	for(int i=0;i<N;i++) sigma[i]=(double*)malloc(sizeof(double)*N);

	// double sigmaInv[N][N];
	double **sigmaInv=(double**)malloc(sizeof(double*)*N);
	for(int i=0;i<N;i++) sigmaInv[i]=(double*)malloc(sizeof(double)*N);

	// double u[M][N];
	double **u=(double**)malloc(sizeof(double*)*M);
	for(int i=0;i<M;i++) u[i]=(double*)malloc(sizeof(double)*N);

	double maxDiff;

	// setting a = d
	// cout<<"A\n";
	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<N;j++){
	// 		// //printf("%d %d \t\t", i,j);
	// 		//printf("%.0f\t", d[i][j] );
	// 	}
	// 	//printf("\n");
	// }
	// //printf("\n");

	// cout<<"transposing A\n";
	#pragma omp parallel num_threads(NUM_THREADS)
	{
		printf("no of threads presently %d\n", omp_get_num_threads());
		#pragma omp for
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				// a[i][j]=dTranspose[i][j];
				e[i][j]=0;
				
			}
			e[i][i]=1;
		}
	}
	// cout<<"A transpose\n";
	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<N;j++){
	// 		// //printf("%d %d \t\t", i,j);
	// 		//printf("%.0f\t", a[i][j] );
	// 	}
	// 	//printf("\n");
	// }

	printf("starting do-while = %f\n", omp_get_wtime()-startTime );

	double time1=0.0;
	double time2=0.0;
	double time3=0.0;
	double time4=0.0;
	double time5=0.0;
	double time6=0.0;
	double time7=0.0;
	double time8=0.0;
	int numIter=0;
	do{
		// for(int i=0;i<N;i++){
		// 	aPrev[i]=a[i][i];
		// }

		#pragma omp parallel num_threads(NUM_THREADS)
		{
			if(omp_get_thread_num()==0) time1 -= omp_get_wtime();
			#pragma omp for
			for(int i=0;i<N;i++){
				// cout<<"i="<<i<<endl;
				for(int j=0;j<N;j++){
					q[i][j]=a[i][j];
					r[i][j]=0;
				}
				// cout<<"init q[i] as ";
				// print_vector(q[i],N);
			}
			if(omp_get_thread_num()==0) time1 += omp_get_wtime();
		}

		for(int i=0;i<N;i++){
			#pragma omp parallel num_threads(NUM_THREADS)
			{
				if(omp_get_thread_num()==0) time2 -= omp_get_wtime();
				#pragma omp for schedule(dynamic,1)
				for(int k=0;k<i;k++){
					double dot=0.0;
					double normSquared=0.0;
					for(int j=0; j<N; j++){
						dot+=a[i][j]*q[k][j];
						normSquared+=q[k][j]*q[k][j];
					}
					if(normSquared<0.0001) continue;
					// for(int j=0; j<N; j++){
					// 	dot+=a[i][j]*q[k][j];
					// 	// normSquared+=q[k][j]*q[k][j];
					// }

					// double dot=dot_product(a[i],q[k],N);
					// double normSquared=dot_product(q[k],q[k],N);
					// cout<<"k="<<k<<"\tdot="<<dot<<"\tnormSq="<<normSquared<<endl;
					double nume = dot/normSquared;
					#pragma omp critical
					{
						for(int j=0;j<N;j++){
							// cout<<"step "<<q[i][j]<<"\t"<<(a[i][j]*dot)/normSquared << "\t"<<q[i][j] - (a[i][j]*dot)/normSquared<<endl;
							q[i][j] -= nume*q[k][j];
						}
					}
				}
				if(omp_get_thread_num()==0) time2 += omp_get_wtime();
			}
		}
		
		// double localSum[NUM_THREADS][N];
		
		// for(int i=0;i<NUM_THREADS;i++){
		// 	for(int j=0;j<N;j++){
		// 		localSum[i][j]=0;
		// 	}
		// }

		// for(int i=0;i<N;i++){
		// 	#pragma omp parallel num_threads(NUM_THREADS)
		// 	{
		// 		#pragma omp for schedule(dynamic,1)
		// 		for(int k=0;k<i;k++){
		// 			double dot=dot_product(a[i],q[k],N);
		// 			double normSquared=dot_product(q[k],q[k],N);
		// 			for(int j=0;j<N;j++){
		// 				localSum[omp_get_thread_num()][j] += (q[k][j]*dot)/normSquared;
		// 			}
		// 		}
		// 	}

		// 	for(int k=0;k<NUM_THREADS && k<i;k++){
		// 		for(int j=0;j<N;j++){
		// 			q[i][j] -= localSum[k][j];
		// 			localSum[k][j]=0;
		// 		}
		// 	}
		// }

		#pragma omp parallel num_threads(NUM_THREADS)
		{
			if(omp_get_thread_num()==0) time3 -= omp_get_wtime();
			#pragma omp for
			for(int i=0;i<N;i++){
				double norm=sqrt(dot_product(q[i],q[i],N));
				if(norm<0.00001) continue;
				for(int j=0;j<N;j++){
					q[i][j]/=norm;
				}
			}
			if(omp_get_thread_num()==0) time3 += omp_get_wtime();
			
			// cout<<"U transpose\n";
			// for(int i=0;i<N;i++){
			// 	for(int j=0;j<N;j++){
			// 		// //printf("%d %d \t\t", i,j);
			// 		//printf("%f\t", q[i][j] );
			// 	}
			// 	//printf("\n");
			// }
			// //printf("\n");

			// cerr<<"normalising\n";
			// for(int i=0;i<N;i++){
			// 	norm=sqrt(dot_product(q[i],q[i],N));
			// 	for(int j=0;j<N;j++){
			// 		q[i][j]/=norm;
			// 	}
			// }
			if(omp_get_thread_num()==0) time4 -= omp_get_wtime();
			#pragma omp for schedule(dynamic,1)
			for(int i=0;i<N;i++){
				for(int j=0;j<=i;j++){
					r[i][j]=dot_product(q[j],a[i],N);
				}
			}
			if(omp_get_thread_num()==0) time4 += omp_get_wtime();

			// #pragma omp for
			// for(int i=0;i<N;i++){
			// 	for(int j=0;j<N;j++){
			// 		if(j>i){
			// 			r[i][j]=0;
			// 			continue;
			// 		} else{
			// 			r[i][j]=dot_product(q[j],a[i],N);
			// 		}
			// 	}
			// }

			// cout<<"R transpose\n";
			// for(int i=0;i<N;i++){
				// for(int j=0;j<N;j++){
					// //printf("%d %d \t\t", i,j);
					//printf("%f\t", r[i][j] );
				// }
				//printf("\n");
			// }
			//printf("\n");
			if(omp_get_thread_num()==0) time5 -= omp_get_wtime();	
			#pragma omp for	
			for(int i=0; i<N; i++) {
				for(int j=0; j<N; j++) {
					aUpdated[i][j]=0; 
					for(int m=0; m<N; m++) {
						aUpdated[i][j]+=q[i][m]*r[m][j];
					}
				}
			}



			// UPDATING E0 INTO E1
			#pragma omp for
			for(int i=0; i<N; i++) {
				for(int j=0; j<N; j++) {
					eUpdated[i][j]=0; 
					for(int m=0; m<N; m++) {
						eUpdated[i][j]+=q[i][m]*e[m][j];
					}
				}
			}
			if(omp_get_thread_num()==0) time5 += omp_get_wtime();


			// if(omp_get_thread_num()==0) time5 -= omp_get_wtime();
			// #pragma omp for
			// for(int i=0;i<N;i++){
			// 	for(int j=0;j<N;j++){
			// 		a[i][j]=0;
			// 		for(int k=0;k<N;k++){
			// 			a[i][j]+=q[i][k]*r[k][j];
			// 		}
			// 	}
			// }
			// for(int i=0;i<N;i++){
			// 	for(int j=0;j<N;j++){
			// 		eUpdated[i][j]=0;
			// 		for(int k=0;k<N;k++){
			// 			eUpdated[i][j]+=q[i][k]*e[k][j];
			// 		}
			// 	}
			// }
			// if(omp_get_thread_num()==0) time5 += omp_get_wtime();
		// }
			if(omp_get_thread_num()==0) time6 -= omp_get_wtime();
			#pragma omp single
			{
				maxDiff=0;
				for(int i=0;i<N;i++){
					if(fabs(a[i][i]-aUpdated[i][i]) > maxDiff)
						maxDiff=fabs(a[i][i]-aUpdated[i][i]);
				}
				for(int i=0;i<N;i++){
					for(int j=0;j<N;j++){
						if(fabs(e[i][j]-eUpdated[i][j]) > maxDiff)
							maxDiff=fabs(e[i][j]-eUpdated[i][j]);
					}
				}
				numIter++;
			}
			if(omp_get_thread_num()==0) time6 += omp_get_wtime();
		// }
					
		// #pragma omp parallel num_threads(NUM_THREADS)
		// {
			if(omp_get_thread_num()==0) time7 -= omp_get_wtime();
			#pragma omp for
			for(int i=0;i<N;i++){
				for(int j=0;j<N;j++){
					e[i][j]=eUpdated[i][j];
					a[i][j]=aUpdated[i][j];
				}
			}
			if(omp_get_thread_num()==0) time7 += omp_get_wtime();
		}

		// cout<<"iter = "<<numIter << "\tmaxDiff="<<maxDiff<<endl;
	}while(maxDiff > 0.000001 && numIter < 6*M);

	printf("inside samarth code\n");

	// printf("\nD0\n\n");
	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<N;j++){
	// 		printf("%f\t", aUpdated[i][j]);
	// 	}
	// 	printf("\n\n");
	// }

	// printf("\nE0\n\n");
	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<N;j++){
	// 		printf("%f\t", eUpdated[i][j]);
	// 	}
	// 	printf("\n\n");
	// }

	// cout<<"numIter = "<<numIter<<endl;
	cout<<"iter = "<<numIter << "\tmaxDiff="<<maxDiff<<endl;

	cout << "time1 " << time1 << endl;
	cout << "time2 " << time2 << endl;
	cout << "time3 " << time3 << endl;
	cout << "time4 " << time4 << endl;
	cout << "time5 " << time5 << endl;
	cout << "time6 " << time6 << endl;
	cout << "time7 " << time7 << endl;

	printf("ending do-while = %f\n", omp_get_wtime()-startTime );

	// cerr<<"transposing\n";
	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<N;j++){
	// 		q2[i][j]=q[j][i];
	// 	}
	// }

/*
	cout<<"Q\n";
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			// printf("%d %d \t\t", i,j);
			printf("%0.f\t", q2[i][j] );
		}
		printf("\n\n");
	}
	// printf("\n");
*/

	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<N;j++){
	// 		if( (i==j && q2[i][j]!=1) || q2[i][j]!=0 )
	// 	}
	// }

// eUpdated has eTranspose now
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			eUpdated[i][j]=e[j][i];
			// printf("%f\t", eUpdated[i][j]);
		}
		// printf("\n");
	}

	// printf("E infinity transpose\n");
	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<N;j++){
	// 		printf("%f\t", e[i][j]);
	// 	}
	// 	printf("\n");
	// }

	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			sigma[i][j]=0;
			sigmaInv[i][j]=0;
		}
		sigma[i][i]=sqrt(a[i][i]);
		if(fabs(sigma[i][i]) > 0.0001 )
			sigmaInv[i][i]=1/sigma[i][i];
	}

	#pragma omp parallel num_threads(NUM_THREADS)
	{
		#pragma omp for
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				u[i][j]=0;
				for(int k=0;k<N;k++){
					u[i][j]+=d[i][k]*eUpdated[k][j];
				}
				u[i][j]*=sigmaInv[j][j];
			}
		}
	}

	// //printf("sigma inv\n");
	// for(int i=0;i<N;i++){
	// 	//printf("%f\t", sigmaInv[i][i]);
	// }
	// //printf("\n");

	// //printf("my u\n");
	// for(int i=0;i<M;i++){
	// 	for(int j=0;j<N;j++){
	// 		//printf("%f\t", u[i][j]);
	// 	}
	// 	//printf("\n");
	// }

	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			(*U)[N*i+j]=eUpdated[i][j];
		}
	}

	for(int i=0;i<N;i++){
		for(int j=0;j<M;j++){
			(*V_T)[M*i+j]=u[j][i];
		}
	}
	for(int i=N;i<M;i++){
		for(int j=0;j<M;j++){
			(*V_T)[M*i+j]=0;
		}
	}
	
	for(int i=0;i<N;i++){
		*(*SIGMA +i)=sigma[i][i];
	}


	// printf("U to be ret\n");
	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<N;j++){
	// 		printf("%f\t", (*U)[N*i+j]);
	// 	}
	// 	printf("\n");
	// }

	// printf("SIGMA\n");
	// for(int i=0;i<N;i++){
	// 	printf("%f\n", *(*SIGMA +i));
	// }

	// printf("V_T to be ret\n");
	// for(int i=0;i<M;i++){
	// 	for(int j=0;j<M;j++){
	// 		printf("%f\t", (*V_T)[M*i+j]);
	// 	}
	// 	printf("\n");
	// }

	// double sam=0;
	// //printf("V_T to be returned\n");
	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<M;j++){
	// 		sam=dTranspose[i][j];
	// 		for(int k=0;k<N;k++){
	// 			sam -= eUpdated[i][k] * sigma[k][k] * (*V_T)[M*k+j];
	// 		}
	// 		// if(fabs(sam) >= 0.0001 )
	// 			//printf("i=%d j=%d temp=%f\n", i,j, sam);
	// 	}
	// 	// //printf("\n");
	// }
	printf("ending svd = %f\n", omp_get_wtime()-startTime );

}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    //printf("called pca\n");

	// printf("SIGMA\n");
 //    for(int i=0;i<N;i++){
 //    	printf("%f\n", SIGMA[i]);
 //    }

    float* sigmaSorted = (float*)malloc(sizeof(float)*N);
    float sum=0, ret=retention/100.0;

    for(int i=0;i<N;i++){
    	sigmaSorted[i]=SIGMA[i]*SIGMA[i];
    	sum+=sigmaSorted[i];
    }

    // //printf("sigmaSorted non normalised\n");
    // for(int i=0;i<N;i++){
    // 	//printf("%f\n", sigmaSorted[i]);
    // }
    // //printf("sum=%f\n", sum);

    // printf("sigmaSorted\n");
    // for(int i=0;i<N;i++){
    // 	printf("%f\n", sigmaSorted[i]);
    // }

    for(int i=0;i<N;i++){
    	sigmaSorted[i]/=sum;
    }

    // printf("sigmaSorted\n");
    // for(int i=0;i<N;i++){
    // 	printf("%f\n", sigmaSorted[i]);
    // }

    sort(sigmaSorted,sigmaSorted + N,greater<float>());

    // printf("ret = %f\n", ret);

	// printf("sigmaSorted\n");
 //    for(int i=0;i<N;i++){
 //    	printf("%f\n", sigmaSorted[i]);
 //    }

    for(int i=0;i<N;i++){
    	ret-=sigmaSorted[i];
    	if(ret<=0){
    		*K=i+1;
    		break;
    	}
    }

    // //printf("sigmaSorted\n");
    // for(int i=0;i<N;i++){
    // 	//printf("%f\n", sigmaSorted[i]);
    // }

    // //printf("SIGMA\n");
    // for(int i=0;i<N;i++){
    // 	//printf("%f\n", (SIGMA[i]*SIGMA[i])/sum);
    // }

    // float w[N][K];
	float **w=(float**)malloc(sizeof(float*)*N);
	for(int i=0;i<N;i++) w[i]=(float*)malloc(sizeof(float)* *K);

	int j;
    for(int i=0;i<*K;i++){
    	for(j=0;j<N;j++){
    		if(fabs(sigmaSorted[i]-(SIGMA[j]*SIGMA[j]/sum) )<0.0001 )
    			break;
    	}
    	for(int k=0;k<N;k++){
    		w[k][i]=U[N*k+j];
    	}
    }

    // printf("\nw\n\n");
    // for(int i=0;i<N;i++){
    // 	for(int j=0;j<*K;j++){
    // 		printf("%f\t", w[i][j]);
    // 	}
    // 	printf("\n");
    // }

    *D_HAT = (float*)malloc(sizeof(float)*M* *K);

    #pragma omp parallel num_threads(NUM_THREADS)
    {
    	#pragma omp for
	    for(int i=0;i<M;i++){
	    	for(int j=0;j<*K;j++){
	    		(*D_HAT)[*K *i + j] =0;
	    		for(int l=0;l<N;l++){
		    		(*D_HAT)[*K *i+j] += D[N*i + l]*w[l][j];
	    		}
	    	}
	    }
	}
    printf("K=%d\n",*K);

    // printf("\nD_HAT\n\n");
    // for(int i=0;i<M;i++){
    // 	for(int j=0;j<*K;j++){
    // 		printf("%.6f\t", (*D_HAT)[*K *i+j]);
    // 	}
    // 	printf("\n");
    // }
}