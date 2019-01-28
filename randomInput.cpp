#include<iostream>
#include<stdlib.h>
#include<time.h>
using namespace std;

int main(){
	srand ( time(NULL) );

	int k,n;
	cout<<"Enter K\n";
	cin>>k;
	
	cout<<"Enter number of points\n";
	cin>>n;

	cout<<k<<endl<<n<<endl;
	for(int i=0;i<n;i++){
		for(int i=0;i<3;i++){
			cout<<rand()%100<<" ";
		}
		cout<<endl;
	}

	return 0;
}