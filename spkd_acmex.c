/*------------------------------------------------------------------------------------------------*/
/* (MEX FILE) Calculates the distance between spike trains using the Jon Victor spike time metric */
/*                                                                                                */
/*  MATLAB USAGE:   d = spkd_acmex(S1,S2,cost)                                                    */
/*------------------------------------------------------------------------------------------------*/


#include "mex.h"
#include <stdio.h>
/* Input Arguments */
//#define S1  prhs[0]
//#define S2  prhs[1]
//#define cost  prhs[2]
/* Output Arguments */
//#define D plhs[0]

#if !defined(MAX)
#define MAX(A,B)  ((A) > (B) ? (A) : (B))
#endif
#if !defined(MIN)
#define MIN(A,B)  ((A) < (B) ? (A) : (B))
#endif

/*-------------------------------------------------------------------------------------*/
/* Numerical Recipes related definitions */
#define NR_END 1
#define FREE_ARG char*
#define EPS .0000001
#define MAXCOST 10000

void nrerror(char error_text[])
/* numerical recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;
  
  v = (double *) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in dvector()");
  return v-nl+NR_END;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
  double **m;
  /* allocate pointers to rows */
  m=(double **)malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in dmatrix()");
  m += NR_END;
  m -= nrl;
  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  /* return pointer to array of pointers to rows */
  return m;
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

double dabs(double x)
{
  if (x<0) return -x;
  else return x;
}

static void getdist(
                    double    *d,
                    unsigned int nspi,
                    double tli[],
                    unsigned int nspj,
                    double tlj[],
                    double  cost
                    )
{
  double **scr;
  unsigned int i,j;
  
  if (cost < EPS) {
    *d = dabs((double)nspi-(double)nspj);
  }
  else if (cost >= MAXCOST)
    *d = (double)(nspi+nspj);
  else {
    /* Initialize margins with cost of adding a spike */
    scr = dmatrix(0,nspi,0,nspj);
    for (i=0; i<=nspi;i++)
      scr[i][0] = i;
    for (j=0; j<=nspj;j++)
      scr[0][j] = j;
    if (nspi!=0 && nspj!=0){
      /* The heart of the algorithm */
      for (i=1;i<=nspi;i++)
        for (j=1;j<=nspj;j++)
          scr[i][j]=MIN(MIN(scr[i-1][j]+1,scr[i][j-1]+1),scr[i-1][j-1]+cost*dabs(tli[i-1]-tlj[j-1]));
    }
    *d = scr[nspi][nspj];
    free_dmatrix(scr,0,nspi,0,nspj);
  }
  return;
}
/*-------------------------------------------------------------------------------------*/
/*                           MATLAB GATEWAY FUNCTION                                   */
/*-------------------------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *tli, *tlj;
  double *d;
  double cost;
  long i1,i2,j1,j2;
  unsigned int nspi, nspj;
  
  if (nrhs != 3){
    mexErrMsgTxt("Exactly 3 input arguments required:  (spiketimes1, spiketimes2, cost)");    
  }
  if (nlhs > 1){
    mexErrMsgTxt("Too many output arguments!");
  }
  
  /* Get spike train times and their lengths (allow spike times to be a row or column) */
  tli = mxGetPr(prhs[0]);   // spike train 1  (vector, pointer)
  tlj = mxGetPr(prhs[1]);   // spike train 2  (vector, pointer)
  cost = mxGetScalar(prhs[2]);  // cost parameter (scalar)
  i1 = mxGetM(prhs[0]); 
  i2 = mxGetN(prhs[0]);
  nspi = (unsigned int)MAX(i1,i2);  // ensures that either a row or column vector can be given as input
  j1 = mxGetM(prhs[1]); 
  j2 = mxGetN(prhs[1]);
  nspj = (unsigned int)MAX(j1,j2);
  
  /* Create a matrix (here a scalar) for the output argument */
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  /* Assign pointers to the output parameters */
  d = mxGetPr(plhs[0]);   
  
  getdist(d, nspi, tli, nspj, tlj, cost);
  return;
}

