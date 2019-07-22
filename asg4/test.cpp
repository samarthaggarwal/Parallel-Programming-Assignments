// #include "lab4_mpi.h"
#include<iostream>
using namespace std;

#include<stdlib.h>
// #include <malloc.h>
// #include "mpi.h"

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

	// cout<<"pi_y="<<pi_y<<endl;
	// cout<<"phi_y\n";
	// for(int i=0;i<pi_y ; i++){
	// 	cout<<phi_y[i]<<" ";
	// }
	// cout<<"\n";
	// cout<<"n="<<n<<" m="<<m<<endl;

	int *pos = np_text_analysis(text, p_dash, phi_y, n, 2*p - 1);
	
	// cout<<"printing pos in np_text\n";
	// int i;
	// for(i=0; pos[i]!=-1 ; i++){
	// 	cout<<pos[i]<<" ";
	// }
	// cout<<pos[i];
	// cout<<endl;

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
	// for(int i=0 ; i < len_u2v ; i++ ) {
	// 	cout<<u2v[i]<< " ";
	// }
	// cout<<endl;

	int *M = (int *)malloc(sizeof(int) * n);
	bool found_i , found_u2v;
	for(int i=0 ; i < n ; i++ ){
		M[i] = 0;
		found_i = false;
		for(int j=0 ; pos[j] != -1 ; j++){
			if(pos[j] == i){
				found_i = true;
				break;
			}
		}

		found_u2v = true;
		for(int j=0 ; j < len_u2v  && found_u2v ; j++ ){
			if(text[i + j] != u2v[j] ){
				found_u2v = false;
				// cout<<"not found for i="<<i<<endl;
			}
		}

		if( found_i && found_u2v){
			M[i] = 1;
		}
	}

	// for(int i=0; i< n ; i++){
	// 	cout<<M[i]<<" ";
	// }
	// cout<<endl;

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
				if(M[temp2] == 1){
					count++;
					if(count>=k-1){
						C[i][j] = 1;
						count=0;
						break;
					}
				} else {
					count=0;
				}
			}
		}
	}

	// for(int i =0 ; i< p ; i++) {
	// 	for(int j=0;j<ceil_np ; j++){
	// 		cout<<C[i][j];
	// 	}
	// 	cout<<endl;
	// }

	// int* match1 = (int*)malloc( sizeof(int) * (n-m+1) );
	// for(int j = 0 ; j < n-m+1 ; j++ ){
	// 	match1[j] = C[j%p][j/p];
	// }

	// for(int i=0;i<n-m+1;i++){
	// 	cout<<match1[i];
	// }
	// cout<<endl;

	int indexCount=0;
	for(int j = 0 ; j < n-m+1 ; j++ ){
		if(C[j%p][j/p] == 1)
			indexCount++;
	}
	*match = (int*)malloc( sizeof(int) * indexCount );
	int index = 0;
	for(int j = 0 ; j < n-m+1 ; j++ ){
		if(C[j%p][j/p]==1){
			// cout<<"j="<<j<<endl;
			// cout<<"index="<<index<<endl;
			(*match)[index] = j;
			index++;
		}
	}

	// cout<<"printing match\n";
	// for(int i=0;i<indexCount ; i++){
	// 	cout<<(*match)[i]<<" ";
	// }
	// cout<<"===="<<endl;

	return indexCount;
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */

/*
To be implemented
Arguments:
    n     : integer containing number of characters (length) of text 
            in which patterns are to be searched (input)
    text  : character array containing text in which patterns are to be searched (input)
    num_patterns: integer containing #patterns to be matched in the text (input)
    m_set : integer array containing length of patterns (input)
            #elements in m_set = num_patterns
            --------------------------------------------------------------------------
            | len(pattern[0]) | len(pattern[1]) | ... | len(pattern[num_patterns-1]) |
            --------------------------------------------------------------------------
    p_set : integer array containing length of period in patterns (input)
            #elements in p_set = num_patterns
            -----------------------------------------------------------------------------------
            | period(pattern[0]) | period(pattern[1]) | ... | period(pattern[num_patterns-1]) |
            -----------------------------------------------------------------------------------
    pattern_set : array of character array containing patterns to be matchede (input)
            #elements in pattern_set = num_patterns
            -----------------------------------------------------------
            | pattern[0] | pattern[1] | ... | pattern[num_patterns-1] |
            -----------------------------------------------------------
    match_counts : 1D integer array containing number of matches of each pattern in text (output)
            #elements in ocuurance_count = num_patterns
            -----------------------------------------------------------------------------------------
            | #matches(pattern[0]) | #matches(pattern[1]) | ... | #matches(pattern[num_patterns-1]) |
            -----------------------------------------------------------------------------------------
    matches : 1D array of integers containing list of all matches (start index of) of pattern_i in text (output)
            consider index of text starting from 0 (not 1)
            #elements in matches = sum(match_counts)
            -------------------------------------------------------------------------------------------------
            | match(pattern[0])[0] | match(pattern[0])[1] | ... | match(pattern[0])[#matches(pattern[0])-1] |
            -------------------------------------------------------------------------------------------------
            | match(pattern[1])[0] | match(pattern[1])[1] | ... | match(pattern[1])[#matches(pattern[1])-1] |
            -------------------------------------------------------------------------------------------------
            | ... ... ... ... ... ... | match(pattern[num_patterns-1])[#matches(pattern[num_patterns-1])-1] |
            -------------------------------------------------------------------------------------------------
*/

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

	// for(int i=0;i<16;i++)
	// 	cout<<text[i];
	// cout<<endl;
	
	// for(int i=0;i<7;i++)
	// 	cout<<*(*pattern_set + i);
	// cout<<endl;

	*match_counts = (int *)malloc( sizeof(int) * num_patterns );
	int **matches_2d = (int **)malloc( sizeof(int*) * num_patterns );
	int *match;
	int num_matches, total_matches=0;
	for(int i=0 ; i < num_patterns ; i++ ){
		num_matches = p_text_analysis(text, pattern_set[i] , p_set[i] , n , m_set[i] , &match);
		// cout<<"printing match in yoyo\n";
		// for(int i=0;i<num_matches ; i++){
		// 	cout<<match[i]<<" ";
		// }
		// cout<<"===="<<endl;

		*(*match_counts + i) = num_matches;
		total_matches += num_matches;
		matches_2d[i] = (int *)malloc( sizeof(int) * num_matches );
		for(int j=0 ; j < num_matches ; j++) {
			matches_2d[i][j] = match[j];
		}
	}

	// cout<<"printing matches 2d\n";
	// for(int i=0; i< num_patterns ; i++){
	// 	for(int j=0;j<num_matches ; j++){
	// 		cout<<matches_2d[i][j]<<" ";
	// 	}
	// 	cout<<endl;
	// }

	*matches = (int *)malloc( sizeof(int) * total_matches );
	int index = 0;
	for(int i=0 ; i < num_patterns ; i++ ){
		for(int j=0 ; j < *(*match_counts + i) ; j++ ){
			*(*matches + index) = matches_2d[i][j];
			index++;
		}
	}

}

int main(){
// 	char text[16] = {'b','a','b','a','a','b','a','b','a','a','b','a'};
// 	char pat[5] = {'a','b','a','a','b'};
// 	int phi[3] = {0,0,1};

// 	int *match = np_text_analysis(text , pat, phi , 12 , 5);
// 	cout<<"====";
// 	for(int i=0 ; match[i] != -1 ; i++){
// 		cout<<match[i] << " ";
// 	}
// 	cout<<endl;
// 	return 0;
// }

	int n=16;
	int num_patterns = 1;

	int *m_set = (int*)malloc( sizeof(int) * num_patterns);
	m_set[0] = 7;

	int *p_set = (int*)malloc( sizeof(int) * num_patterns);
	p_set[0] = 2;
	
	char text[16] = {'b','a','b','a','b','a','b','a','b','a','b','a','a','b','a','b'};
	char **pattern_set = (char**)malloc( sizeof(char*) * num_patterns);
	*pattern_set = (char *)malloc( sizeof(char) * 7);
	char pattern[7] = {'a','b','a','b','a','b','a'};

	for(int i=0;i<7;i++)
		*(*pattern_set + i) = pattern[i];
	// cout<<endl;

	// for(int i=0;i<16;i++)
		// cout<<text[i];
	// cout<<endl;
	
	// for(int i=0;i<7;i++)
		// cout<<pattern[i];
	// cout<<endl;

	// for(int i=0;i<7;i++)
		// cout<<*(*pattern_set + i);
	// cout<<endl;

	int *match_counts;
	int *matches;

	// cout<<"calling per_pat_match\n";
	periodic_pattern_matching( n, text, num_patterns, m_set, p_set, pattern_set, &match_counts, &matches);

	int index=0;
	cout<<"match_counts, matches\n";
	for(int i=0 ; i < num_patterns ; i++){
		cout<<match_counts[i]<<"\n";
		for(int j=0 ; j < match_counts[i] ; j++){
			cout<<"index="<<index<<" value="<<matches[index]<<"\n";
			index++;
		}
	}

	return 0;
}
