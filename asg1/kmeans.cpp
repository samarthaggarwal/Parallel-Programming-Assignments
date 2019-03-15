#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<climits>
#include<algorithm>
using namespace std;

int assignPoints(int n, int k, int points[][4], float means[][4]){ // assigns points to means as 4th dimension
	
	float temp, tempDist, minDist, minDistIndex;
	int numChanges = 0;

	for(int i=0;i<n;i++){
		minDist=INT_MAX;
		minDistIndex=-1;
		for(int j=0;j<k;j++){

			tempDist=0;
			for(int k=0;k<3;k++){
				temp=points[i][k]-means[j][k];
				tempDist += temp*temp;
			}

			if(tempDist<minDist){
				minDist=tempDist;
				minDistIndex=j;
			}
		}

		if(minDistIndex!=points[i][3]){
			points[i][3]=minDistIndex;
			numChanges++;
		}
	}

	return numChanges;
}

void recomputeMeans(int n, int k, int points[][4], float means[][4]){ // recompute means for each cluster
	int meanIndex;

	for(int i=0;i<k;i++){
		means[i][3]=0;
	}

	for(int i=0;i<n;i++){
		meanIndex = points[i][3];
		for(int j=0;j<3;j++){
			means[meanIndex][j] = 0;
		}
	}

	for(int i=0;i<n;i++){
		meanIndex = points[i][3];
		for(int j=0;j<3;j++){
			means[meanIndex][j]+=points[i][j];
		}
		means[meanIndex][3]++;
	}

	for(int i=0;i<k;i++){
		if(means[i][3]==0)
			continue;
		for(int j=0;j<3;j++){
			means[i][j]/=means[i][3];
		}
	}

	return;
}

void printPoints(int points[][4], int n){
	cout<<"\nPoints\n";
	for(int i=0;i<n;i++){
		for(int j=0;j<4;j++){
			cout<<points[i][j]<<" ";
		}
		cout<<endl;
	}
}

void printMeans(float means[][4], int k){
	cout<<"\nK means\n";
	for(int i=0;i<k;i++){
		for(int j=0;j<4;j++){
			cout<<means[i][j]<<" ";
		}
		cout<<endl;
	}
}

int main(int argc, char *argv[]){
	srand (time(NULL));

	int k,n, maxIterations, thresNumChanges;
	cout<<"Enter K\n";
	cin>>k;
	
	cout<<"Enter number of points\n";
	cin>>n;
	int points[n][4];// 4th dim - cluster number
	float means[k][4];// 4th dim - no. of points in cluster

	//reading points
	for(int i=0;i<n;i++){
		for(int j=0;j<3;j++){
			cin>>points[i][j];
		}
	}

	random_shuffle(&points[0],&points[n]);
	// initialising means
	for(int i=0; i<k; i++){
		for(int j=0;j<3;j++){
			// means[i][j]=rand()%50;
			means[i][j]=points[i][j];
		}
	}

	maxIterations = 200;
	thresNumChanges = 0;
	for(int i=0;i<maxIterations;i++){
		// cout<<"\n\niter "<<i+1<<endl;
		// printPoints(points,n);
		// printMeans(means,k);
	
		if(assignPoints(n,k,points,means) <= thresNumChanges){
			cout<<"numIter = "<<i+1<<endl;
			break;
		}
		recomputeMeans(n,k,points,means);
	}

	// printPoints(points,n);
	printMeans(means,k);

	return 0;
}