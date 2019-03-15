#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<climits>
#include<algorithm>
#include<pthread.h>
#include<omp.h>
using namespace std;

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

void assign_points_t(){ // assigns points to means as 4th dimension

	point *points = pointsPtr;
	mean *means = meansPtr;

	#pragma omp parallel for
	for(int i=0 ;i<n ;i++){

		float minDist=INT_MAX;
		float minDistIndex=-1;
		float temp, tempDist;
		for(int j=0;j<k;j++){

			// cerr<<"point"<<i+1<<"\t"<<points[i].x<<" "<<points[i].y<<" "<<points[i].z<<"\t\tmean"<<j+1<<" "<<means[i].x<<" "<<means[i].y<<" "<<means[i].z<<"\t";


			tempDist=0;
			temp = points[i].x - means[j].x;
			tempDist += temp*temp;
	 		temp = points[i].y - means[j].y;
			tempDist += temp*temp;
			temp = points[i].z - means[j].z;
			tempDist += temp*temp;

			if(tempDist<minDist){
				minDist=tempDist;
				minDistIndex=j;
			}
		}

		if(minDistIndex!=points[i].cluster){
			points[i].cluster=minDistIndex;
			// numChanges++;
		}
	}
	// for(int i=0;i<10;i++){
	// 	cerr<<"threadID="<<omp_get_thread_num()<<"\ti="<<i;
	// 	cerr<<"hello\n";
	// }
	#pragma omp barrier
	// cerr<<endl;
	// }

	// return numChanges;
}

void recompute_means(){ // recompute means for each cluster
	
	point *points = pointsPtr;
	mean *means = meansPtr;
	int meanIndex;

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
	// return;
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

int main(int argc, char *argv[]){
	// srand (time(NULL));
	srand(3);

	int maxIterations, thresNumChanges;
	cout<<"Enter K\n";
	cin>>k;
	
	cout<<"Enter number of points\n";
	cin>>n;
	point points[n];// 4th dim - cluster number
	mean means[k];// 4th dim - no. of points in cluster
	
	numThreads = 8;
	lengthPerThread = (n / numThreads) + 1;
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(numThreads);

	pointsPtr = points;
	meansPtr = means;

	//reading points
	for(int i=0;i<n;i++){
		cin>>points[i].x;
		cin>>points[i].y;
		cin>>points[i].z;
		points[i].cluster=0;
	}

	random_shuffle(&points[0],&points[n]);
	// initialising means
	for(int i=0; i<k; i++){
		// means[i][j]=rand()%50;
		means[i].x=points[i].x;
		means[i].y=points[i].y;
		means[i].z=points[i].z;
	}

	double start, end;
	start = omp_get_wtime();

	// pthread_t threads[numThreads];
	// pthread_mutex_init(&lock, NULL);

	maxIterations = 100;
	thresNumChanges = 0;
	for(int i=0;i<maxIterations;i++){
		// cout<<"iter "<<i+1<<endl;
		// printPoints(points);
		// printMeans(means);

		assign_points_t();
		// printPoints();

		// if(assign_points_t(points,means) <= thresNumChanges){
		// 	cout<<"ended at numIter = "<<i+1<<endl;
		// 	break;
		// }
		recompute_means();
	}

	end = omp_get_wtime();
	// printPoints();
	printMeans();
	cout<<"time = "<<end-start<<endl;

	return 0;
}