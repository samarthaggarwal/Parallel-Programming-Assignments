#include<iostream>
#include<stdlib.h>
#include<time.h>
using namespace std;

int main(int argc, char *argv[]){
	srand (time(NULL));

	int k,n;
	cout<<"Enter K\n";
	cin>>k;
	
	cout<<"Enter number of points\n";
	cin>>n;
	int points[n][3];
	float means[k][3];

	for(int i=0;i<n;i++){
		for(int j=0;j<3;j++){
			cin>>points[i][j];
		}
	}

	for(int i=0; i<k; i++){
		for(int j=0;j<3;j++){
			means[i][j]=rand()%100;
		}
	}

	cout<<"points\n";
	for(int i=0;i<n;i++){
		for(int j=0;j<3;j++){
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