# include <R.h>
void hommola(
        int *n,      
        int *m,       
        int *N,         
        double *ht,     
        double *pt,    	         
	int *posh,
        int *posp,	         
	double *dh,
        double *dp
){         

    int i, j, cnt;
  
    cnt=0;
    for(i=0; i<N[0]; i++){
      for(j=(i+1); j<N[0]; j++){
        dh[cnt] = ht[posh[i]*n[0]+posh[j]];
        dp[cnt] = pt[posp[i]*m[0]+posp[j]];
        cnt ++;
      }
    }

}

