#include<iostream>
#include <pthread.h>
#include <stdlib.h>
using namespace std;

#define N 5

void *worker_thread(void *arg)
{
    printf("This is worker_thread() 1 #%ld\n", (long)arg);
    printf("This is worker_thread() 2 #%ld\n", (long)arg);
    printf("This is worker_thread() 3 #%ld\n", (long)arg);
    printf("This is worker_thread() 4 #%ld\n", (long)arg);
    printf("This is worker_thread() 5 #%ld\n", (long)arg);
    printf("This is worker_thread() 6 #%ld\n", (long)arg);
    pthread_exit(NULL);
}

int main()
{
    pthread_t my_thread[N];
    int ret;

    cout<<"In main: creating thread\n";

    for(long id=0;id<=N;id++){
	    ret =  pthread_create(&my_thread[id], NULL, &worker_thread, (void *)id);
	    if(ret != 0) {
	        cout<<"Error: pthread_create() failed\n";
	        exit(EXIT_FAILURE);
	    }
	}
	
    pthread_exit(NULL);
}