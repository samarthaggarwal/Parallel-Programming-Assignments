#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<climits>
#include<algorithm>
#include<pthread.h>
#include<omp.h>
#include<cmath>
using namespace std;

// pthread_mutex_t lock;
int n;
int k;
int lengthPerThread;
int numThreads;

struct point{
	int x,y,z;
	int cluster;
};

struct mean{
	int count;
	float x,y,z;
};

point *pointsPtr;
mean *meansPtr;

void assign_points_t(int *tid){ // assigns points to means as 4th dimension

	point *points = pointsPtr;
	mean *means = meansPtr;

	float temp, tempDist, minDist, minDistIndex;
	// int numChanges = 0;
	// int *tt = (int*) tid;
	// int t = *tt;
	int t = *tid;

	for(int i=t*lengthPerThread ; i<(t+1)*lengthPerThread && i<n ; i++){
		minDist=INT_MAX;
		minDistIndex=-1;
		for(int j=0;j<k;j++){

			// cerr<<"point"<<i+1<<"\t"<<points[i].x<<" "<<points[i].y<<" "<<points[i].z<<"\t\tmean"<<j+1<<" "<<means[i].x<<" "<<means[i].y<<" "<<means[i].z<<"\t";


			tempDist=0;
			temp = points[i].x - means[j].x;
			tempDist += temp*temp;
	 		temp = points[i].y - means[j].y;
			tempDist += temp*temp;
			temp = points[i].z - means[j].z;
			tempDist += temp*temp;

			// cerr << tempDist <<endl;
			
			// for(int k=0;k<3;k++){
			// 	temp=points[i][k]-means[j][k];
			// 	tempDist += temp*temp;
			// }

			if(tempDist<minDist){
				minDist=tempDist;
				minDistIndex=j;
			}
		}

		if(minDistIndex!=points[i].cluster){
			points[i].cluster=minDistIndex;
			// numChanges++;
		}

		// cerr << minDistIndex<< " ";
	}
	// cerr<<endl;

	// return NULL;
}

void recompute_means(){ // recompute means for each cluster
	
	point *points = pointsPtr;
	mean *means = meansPtr;
	int meanIndex;
	// int *tt = (int*) tid;
	// int t = *tt;

	for(int i=0;i<k;i++){
		means[i].count=0;
	}

	for(int i=0;i<n;i++){
		meanIndex = points[i].cluster;
		means[meanIndex].x=0;
		means[meanIndex].y=0;
		means[meanIndex].z=0;
		// for(int j=0;j<3;j++){
		// 	means[meanIndex][j] = 0;
		// }
	}

	for(int i=0;i<n;i++){
		meanIndex = points[i].cluster;
		means[meanIndex].x += points[i].x;
		means[meanIndex].y += points[i].y;
		means[meanIndex].z += points[i].z;
		
		// for(int j=0;j<3;j++){
		// 	means[meanIndex][j]+=points[i][j];
		// }
		means[meanIndex].count++;
	}

	for(int i=0;i<k;i++){
		if(means[i].count==0)
			continue;
		means[i].x/=means[i].count;
		means[i].y/=means[i].count;
		means[i].z/=means[i].count;
		
		// for(int j=0;j<3;j++){
		// 	means[i][j]/=means[i][3];
		// }
	}

	return;
}

void printPoints(){
	point *points = pointsPtr;

	cout<<"\nPoints\n";
	for(int i=0;i<n;i++){
		cout<<points[i].x<<" "<<points[i].y<<" "<<points[i].z<<" "<<points[i].cluster<<endl;
		// for(int j=0;j<4;j++){
		// 	cout<<points[i][j]<<" ";
		// }
		// cout<<endl;
	}
}

void printMeans(){
	mean *means = meansPtr;

	cout<<"\nK means\n";
	for(int i=0;i<k;i++){
		cout<<means[i].x<<" "<<means[i].y<<" "<<means[i].z<<" "<<means[i].count<<endl;
		// for(int j=0;j<4;j++){
		// 	cout<<means[i][j]<<" ";
		// }
		// cout<<endl;
	}
}

double cost(){
	int meanIndex;
	float temp , tempDist;
	double cost = 0;

	for(int i=0;i<n;i++){
		tempDist = 0;
		meanIndex = pointsPtr[i].cluster;
		temp = pointsPtr[i].x - meansPtr[meanIndex].x;
		tempDist += temp * temp;
		temp = pointsPtr[i].y - meansPtr[meanIndex].y;
		tempDist += temp * temp;
		temp = pointsPtr[i].z - meansPtr[meanIndex].z;
		tempDist += temp * temp;
		cost += tempDist;
	}

	return sqrt(cost);
}

void kmeans_omp(int num_threads, int N, int K, int* data_points, int** data_point_cluster, float** centroids, int* num_iterations){

	// double start, end;
	// start = omp_get_wtime();
	
	// srand (time(NULL));
	// srand(2);

	int maxIterations = 100;
	// cout<<"Enter K\n";
	// cin>>k;
	k=K;
	
	// cout<<"Enter number of points\n";
	// cin>>n;
	n=N;

	numThreads = num_threads;
	lengthPerThread = (n/numThreads) + 1;

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(numThreads);

	*centroids = (float *)malloc(sizeof(float) * maxIterations * 3 * k);
	// float* outputMeans = *centroids;
	int base = 0;

	// point points[n];
	// mean means[k];

	point* points = (point *) malloc (sizeof(point) * n);
	mean* means = (mean *) malloc (sizeof(mean) * k);

	pointsPtr = points;
	meansPtr = means;

	//reading points
	for(int i=0;i<n;i++){
		points[i].x = *(data_points + 3*i);
		points[i].y = *(data_points + 3*i + 1);
		points[i].z = *(data_points + 3*i + 2);
	}
	
	// for(int i=0;i<n;i++){
	// 	cin>>points[i].x;
	// 	cin>>points[i].y;
	// 	cin>>points[i].z;
	// }

	// random_shuffle(&points[0],&points[n]);
	// initialising means
	for(int i=0; i<k; i++){
		// means[i][j]=rand()%50;
		means[i].x=points[i].x;
		means[i].y=points[i].y;
		means[i].z=points[i].z;
	}


	// pthread_t threads[numThreads];
	// pthread_mutex_init(&lock, NULL);
	// int tid[numThreads];
	int *tid = (int *) malloc (sizeof (int) * numThreads);
	for (int t = 0; t < numThreads; ++t){
		tid[t]=t;
	}

	// thresNumChanges = 0;
	*num_iterations = 0;
	for(int iter = 0; iter < maxIterations ; iter++ ){

		*num_iterations = *num_iterations + 1;

		for(int i=0;i<k;i++){
			*(base + (*centroids) + 3*i + 0) = (float)means[i].x;
			*(base + (*centroids) + 3*i + 1) = (float)means[i].y;
			*(base + (*centroids) + 3*i + 2) = (float)means[i].z;
		}
		base += 3*k;

		#pragma omp parallel for
		for(int t=0;t<numThreads;t++){
			// cerr << t;
			// cerr << omp_get_thread_num() <<"\t";
			// pthread_create(&threads[t], NULL, assign_points_t, &tid[t]);
			assign_points_t(&tid[t]);
		}
		// cerr<<endl;

		#pragma omp barrier
		// for (int t=0; t<numThreads; t++){
  // 			pthread_join(threads[t], NULL);
		// }

		// assign_points_t();
		// printPoints();

		// if(assign_points_t(points,means) <= thresNumChanges){
		// 	cout<<"ended at numIter = "<<i+1<<endl;
		// 	break;
		// }
		mean temp_means[k];
		for(int i=0;i<k;i++){
			temp_means[i] = means[i];
		}

		recompute_means();
		
		int i;
		for(i=0;i<k;i++){
			if( fabs(temp_means[i].x - means[i].x) > 0.001 || fabs(temp_means[i].y - means[i].y) > 0.001 || fabs(temp_means[i].z - means[i].z) > 0.001 ){
				break;
			}
		}

		if(i>=k){
			break;
		}
	}

	// printMeans();
	// printPoints();

	// cout << *num_iterations<<endl;
	// writing centroid after last recomputation
	for(int i=0;i<k;i++){
		*(base + (*centroids) + 3*i + 0) = (float)means[i].x;
		*(base + (*centroids) + 3*i + 1) = (float)means[i].y;
		*(base + (*centroids) + 3*i + 2) = (float)means[i].z;
	}
	base += 3*k;

	// writing final points with their cluster IDs
	*data_point_cluster = (int *)malloc(sizeof(int) * 4 * n);
	int* outputPoints = *data_point_cluster;
	for(int i=0;i<n;i++){
		*(outputPoints + 4*i + 0) = points[i].x;
		*(outputPoints + 4*i + 1) = points[i].y;
		*(outputPoints + 4*i + 2) = points[i].z;
		*(outputPoints + 4*i + 3) = points[i].cluster;
	}

	// end = omp_get_wtime();
	// cout<<"time = "<<end-start<<endl;
	// cout<<"cost = "<<cost()<<endl;

	return;
}