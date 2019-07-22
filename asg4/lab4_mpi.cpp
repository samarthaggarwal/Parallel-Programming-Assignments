#include "lab4_mpi.h"
// #include<iostream>
// using namespace std;

#include<stdlib.h>
#include <malloc.h>
#include "mpi.h"

int duel(char *z, char *y, int *phi, int i, int j, int n){
	int k=phi[j-i];
	if(j+k>=n || z[j+k] != y[k])
		return i;
	else
		return j;
}

int *np_text_analysis(char *t, char *p, int *phi, int n, int m){
	int ceilM;
	if(m%2==0){
		ceilM = m/2;
	}else{
		ceilM = (m/2) + 1 ;
	}

	// numT = b+1
	int b;
	if(n%ceilM == 0 )
		b = n/ceilM;
	else
		b = (n/ceilM) + 1;

	// cout<<"n = "<<n<<endl;
	// cout<<"ceilM = "<<ceilM<<endl;
	// cout<<"b = "<<b<<endl;
	int i;
	int potential_positions[b];
	for(int b_i=0 ; b_i < b ; b_i++){

		i = b_i*ceilM;
		for(int j=i+1;j<(b_i+1)*ceilM; j++){
			i = duel(t, p, phi, i, j, n);
		}
		potential_positions[b_i] = i;
	}

	// cout<<"potential positions\n";
	// for(int i=0;i<b;i++)
	// 	cout<<potential_positions[i]<<" ";
	// cout<<endl;

	// 1 extra index to store -1
	int *match_positions = (int *)malloc(sizeof(int) * (b+1) );
	int marker = 0, j;
	for(int b_i = 0 ; b_i < b ; b_i++){
		for(j=0;j<m;j++){
			// cout<<"pot pos ="<<potential_positions[b_i]<< " " << (t[j+ potential_positions[b_i]] != p[j]) <<endl;
			if(t[j+ potential_positions[b_i]] != p[j]){
				// cout<<"i="<<i<<" j="<<j<<" m="<<m<<endl;
				break;
			}
		}
		// cout<<"j="<<j<<" m="<<m<<endl;
		if(j==m){
			// cout<<"adding in array\n";
			match_positions[marker] = potential_positions[b_i];
			marker++;
		}
	}
	// cout<<"marker="<<marker<<endl;
	match_positions[marker] = -1;

	// cout<<"printing match_positions\n";
	// for(int i=0; i<b+1 ; i++){
	// 	cout<<match_positions[i] << " ";
	// }
	// cout<<endl;

	return match_positions;
}

int p_text_analysis(char *text ,char *pattern, int p, int n, int m, int **match){
	
	char p_dash[2*p -1];
	for(int i=0; i<2*p - 1 ; i++){
		p_dash[i] = pattern[i];
	}

	int ceilM;
	if(m%2==0){
		ceilM = m/2;
	}else{
		ceilM = (m/2) + 1 ;
	}

	// witness array for p_dash
	int ceil_p_dash = (2*p-1)/2 + (2*p - 1)%2;
	int pi_y = p<ceil_p_dash ? p : ceil_p_dash;
	int temp;
	int phi_y[pi_y];
	phi_y[0] = 0;
	for(int i=1 ; i < pi_y ; i++ ){
		temp=0;
		while(temp < p && p_dash[i+temp]==p_dash[temp])
			temp++;
		
		phi_y[i]=temp;
	}

	int *pos = np_text_analysis(text, p_dash, phi_y, n, 2*p - 1);
	int k = m/p;

	int len_u2v = (2*p) + m - (k*p);
	char u2v[len_u2v];
	for(int i=0 ; i<p ; i++ ){
		u2v[i] = pattern[i];
		u2v[i+p] = pattern[i];
	}
	for(int i=2*p ; i < len_u2v ; i++ ){
		u2v[i] = pattern[(k-2)*p + i];
	}
	

	int *M = (int *)malloc(sizeof(int) * n);
	for(int i=0; i<n; i++){
		M[i] = 0;
	}
	bool found_i , found_u2v;
	int i;
	for(int iter=0 ; pos[iter]!=-1 ; iter++ ){
		i = pos[iter];
		found_i = true;
		found_u2v = true;
		for(int j=0 ; j < len_u2v  && found_u2v ; j++ ){
			if(text[i + j] != u2v[j] ){
				found_u2v = false;
			}
		}

		if( found_i && found_u2v){
			M[i] = 1;
		}
	}

	int ceil_np;
	if(n%p==0){
		ceil_np = n/p;
	}else{
		ceil_np = (n/p) + 1 ;
	}

	int **C = (int**)malloc(sizeof(int*) * p);
	for(int i=0;i<p;i++){
		C[i] = (int*)malloc(sizeof(int) * (ceil_np+1) ) ;
	}
	int temp1, temp2, count;
	for(int i=0 ; i < p ; i++ ){
		temp1 = i;
		for(int j=0 ; temp1 < n ; j++, temp1 += p ){
			C[i][j] = 0;
			temp2 = temp1;
			count = 0;
			for(int l=0 ; temp2 < n ; l++, temp2 += p){
				if(M[temp2] != 1){
					count=0;
					break;
				} else {
					count++;
					if(count>=k-1){
						C[i][j] = 1;
						count=0;
						break;
					}
				}
			}
		}
	}

	int indexCount=0;
	for(int j = 0 ; j < n-m+1 ; j++ ){
		if(C[j%p][j/p] == 1){
			indexCount++;
		}
	}

	// printf("M=====\n");
	// for(int i=0;i<n;i++){
	// 	printf("%d",M[i]);
	// }
	// cout<<"\n";

	*match = (int*)malloc( sizeof(int) * indexCount );
	int index = 0;
	for(int j = 0 ; j < n-m+1 ; j++ ){
		if(C[j%p][j/p]==1){
			(*match)[index] = j;
			index++;
		}
	}

	for(int i=0; i<p; i++){
		free(C[i]);
	}
	free(pos);
	free(M);
	free(C);

	return indexCount;
}


void periodic_pattern_matching (
		int n, 
		char *text, 
		int num_patterns, 
		int *m_set, 
		int *p_set, 
		char **pattern_set, 
		int **match_counts, 
		int **matches)
{

	int num_processes;
	int process_id;
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_id);

	int patterns_per_proc = num_patterns / num_processes;
	int start, end;
	if(process_id < (num_patterns % num_processes) ){
		start = process_id * (patterns_per_proc + 1);
		end = (process_id+1) * (patterns_per_proc + 1);
	} else {
		int temp = num_patterns % num_processes;
		start = temp * (patterns_per_proc + 1) + (process_id-temp)*patterns_per_proc ;
		end = start + patterns_per_proc;
	}

	int *local_match_counts = (int *)malloc( sizeof(int) * (end - start) );
	int **matches_2d = (int **)malloc( sizeof(int*) * (end - start) );
	int *match;
	int num_matches, total_matches=0;
	for(int i=start ; i < end ; i++ ){
		num_matches = p_text_analysis(text, pattern_set[i] , p_set[i] , n , m_set[i] , &match);
		*(local_match_counts + i - start) = num_matches;
		total_matches += num_matches;
		matches_2d[i-start] = (int *)malloc( sizeof(int) * num_matches );
		for(int j=0 ; j < num_matches ; j++) {
			matches_2d[i-start][j] = match[j];
		}
		free(match);
	}

	int *local_matches = (int *)malloc( sizeof(int) * total_matches );
	int index = 0;
	for(int i=0 ; i < (end - start) ; i++ ){
		for(int j=0 ; j < *(local_match_counts + i) ; j++ ){
			*(local_matches + index) = matches_2d[i][j];
			index++;
		}
	}

	if( process_id == 0){
		*match_counts = (int *)malloc( sizeof(int) * num_patterns );
	}
	int sizes1[num_processes];
	int displ1[num_processes];

	if(process_id == 0){
		for(int i=0;i<num_processes ; i++){
			if( i < (num_patterns % num_processes) ){
				sizes1[i] =  patterns_per_proc + 1;
			} else {
				sizes1[i]  = patterns_per_proc;
			}
		}
		displ1[0] = 0;
		for(int i=1;i<num_processes ; i++ ){
			displ1[i] = displ1[i-1] + sizes1[i-1];
		}
	}

	MPI_Gatherv( local_match_counts, end-start , MPI_INT , *match_counts , sizes1 , displ1 , MPI_INT, 0 , MPI_COMM_WORLD );

	int sizes[num_processes];
	int displ[num_processes];

	if(process_id == 0){
		index=0;
		for(int i=0;i<num_processes ; i++){
			sizes[i]=0;
			for(int j = 0 ; j < sizes1[i] ; j++ ){
				sizes[i] += (*match_counts)[index];
				index++;
			}
		}

		displ[0] = 0;
		for(int i=1;i<num_processes ; i++ ){
			displ[i] = displ[i-1] + sizes[i-1];
		}
	}

	if(process_id == 0){
		*matches = (int *)malloc( sizeof(int) * (displ[num_processes -1] + sizes[num_processes-1]) );
	}
	MPI_Gatherv( local_matches, total_matches , MPI_INT, *matches, sizes, displ, MPI_INT, 0, MPI_COMM_WORLD);
 
}
