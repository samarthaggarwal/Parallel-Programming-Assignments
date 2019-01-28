#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<climits>
using namespace std;

void assign(int n, int k, int points[][4], float means[][3]){ // assigns means to points as 4th dimension
	
	float temp, tempDist, minDist=INT_MAX, minDistIndex=-1;
	
	for(int i=0;i<n;i++){
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

		points[i][3]=minDistIndex;
	}

	return;
}

int main(int argc, char *argv[]){
	srand (time(NULL));

	int k,n;
	cout<<"Enter K\n";
	cin>>k;
	
	cout<<"Enter number of points\n";
	cin>>n;
	int points[n][4];
	float means[k][3];

	for(int i=0;i<n;i++){
		for(int j=0;j<3;j++){
			cin>>points[i][j];
		}
	}

	for(int i=0; i<k; i++){
		for(int j=0;j<3;j++){
			means[i][j]=rand()%50;
		}
	}

	assign(n,k,points,means);

	cout<<"points\n";
	for(int i=0;i<n;i++){
		for(int j=0;j<4;j++){
			cout<<points[i][j]<<" ";
		}
		cout<<endl;
	}

	cout<<"K means\n";
	for(int i=0;i<k;i++){
		for(int j=0;j<3;j++){
			cout<<means[i][j]<<" ";
		}
		cout<<endl;
	}


	return 0;
}