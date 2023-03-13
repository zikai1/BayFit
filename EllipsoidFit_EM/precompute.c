/*
Andriy Myronenko
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mex.h"
#define	max(A, B)	((A) > (B) ? (A) : (B))
#define	min(A, B)	((A) < (B) ? (A) : (B))

void cpd_comp(
		double* x,
		double* y, 
        double* sigma2,
		double* outlier,
        double* P1,
        double* Pt1,
        double* Px,
        double* PM,
    	double* E,
        double* V,
        int N,
		int M,
        int D
        )

{
  int		n, m, d;
  double	ksig, diff, razn, outlier_tmp, outlier_tmp2, sp;
  double	*P, *temp_x;
  
  P = (double*) calloc(M, sizeof(double));
  temp_x = (double*) calloc(D, sizeof(double));
  
  ksig = -2.0 * *sigma2;
  //outlier_tmp=(*outlier*M*pow (-ksig*3.14159265358979,0.5*D))/((1-*outlier)*N); 
  outlier_tmp=(*outlier*M*pow (-ksig*3.14159265358979,0.5*D))/((1-*outlier)**V); 
  //outlier_tmp2=(N*(1-*outlier))/(*outlier*M);
  outlier_tmp2=(*V*(1-*outlier))/(*outlier*M);
  /* printf ("ksig = %lf\n", *sigma2);*/
  /* outlier_tmp=*outlier*N/(1- *outlier)/M*(-ksig*3.14159265358979); */
  
  
  for (n=0; n < N; n++) {
      
      sp=0;
      for (m=0; m < M; m++) {
          razn=0;
          for (d=0; d < D; d++) {
             diff=*(x+n+d*N)-*(y+m+d*M);  
             diff=diff*diff;
             razn+=diff;
          }
          
          *(P+m)=exp(razn/ksig);
          sp+=*(P+m);
      }
      
      *(PM+n)=1/(1+outlier_tmp2*sp);
      sp+=outlier_tmp;
      *(Pt1+n)=1-outlier_tmp/sp;// outlier_free=sp-outlier_tmp; Pt1=outlier_free/sp;
      
      for (d=0; d < D; d++) {
       *(temp_x+d)=*(x+n+d*N)/sp;
      }
         
      for (m=0; m < M; m++) {
         
          *(P1+m)+=*(P+m)/sp;
          
          for (d=0; d < D; d++) {
          *(Px+m+d*M)+= *(temp_x+d)**(P+m);
          }
          
      }
      
   *E +=  -log(sp);     
  }
  *E +=D*N*log(*sigma2)/2-N*log(1-*outlier);
    
  
  free((void*)P);
  free((void*)temp_x);

  return;
}

/* Input arguments */
#define IN_x		prhs[0]
#define IN_y		prhs[1]
#define IN_sigma2	prhs[2]
#define IN_outlier	prhs[3]
#define IN_vol      prhs[4]


/* Output arguments */
#define OUT_P1		plhs[0]
#define OUT_Pt1		plhs[1]
#define OUT_Px		plhs[2]
#define OUT_PM      plhs[3]
#define OUT_E		plhs[4]



/* Gateway routine */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *x, *y, *sigma2, *outlier, *P1, *Pt1, *Px, *PM, *E, *V;
  int     N, M, D;
  
  /* Get the sizes of each input argument */
  N = mxGetM(IN_x);// mxGetM得到行数
  M = mxGetM(IN_y);
  D = mxGetN(IN_x);// mxGetN得到列数
  //V=mxGetScalar(IN_vol);
  
  
  /* Create the new arrays and set the output pointers to them */
  OUT_P1     = mxCreateDoubleMatrix(M, 1, mxREAL);
  OUT_Pt1    = mxCreateDoubleMatrix(N, 1, mxREAL);
  OUT_Px    = mxCreateDoubleMatrix(M, D, mxREAL);
  OUT_E       = mxCreateDoubleMatrix(1, 1, mxREAL);
  OUT_PM     =mxCreateDoubleMatrix(N,1,mxREAL);

    /* Assign pointers to the input arguments */
  x      = mxGetPr(IN_x);
  y       = mxGetPr(IN_y);
  sigma2       = mxGetPr(IN_sigma2);
  outlier    = mxGetPr(IN_outlier);
  V          =mxGetPr(IN_vol);

 
  
  /* Assign pointers to the output arguments */
  P1      = mxGetPr(OUT_P1);
  Pt1      = mxGetPr(OUT_Pt1);
  Px      = mxGetPr(OUT_Px);
  E     = mxGetPr(OUT_E);
  PM    =mxGetPr(OUT_PM);
   
  /* Do the actual computations in a subroutine */
  cpd_comp(x, y, sigma2, outlier, P1, Pt1, Px, PM, E, V, N, M, D);
  
  return;
}


