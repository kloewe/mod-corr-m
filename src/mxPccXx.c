/*----------------------------------------------------------------------------
  File    : mxPccXx.c
  Contents: MEX wrapper for the C-based functions defined in pcc.h
  Author  : Kristian Loewe
----------------------------------------------------------------------------*/
#include "mex.h"
#include "pcc.h"
#include "mxCorrAlloc.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{                                /* --- gateway routine */
  int  N,                        /* number of voxels */
       T,                        /* number of points in time */
       var;                      /* implementation variant */
  int  tilesize = 0;             /* tile size */
  int  nthreads = 0;             /* number of threads */
  REAL *data;                    /* input data array of voxel time series */
  REAL *res;                     /* result */
  int  rc = 1;                   /* return condition */
  
  N    = (int)mxGetN(prhs[0]);   /* retrieve input arguments from MATLAB */
  T    = (int)mxGetM(prhs[0]);
  data = (REAL *)mxGetData(prhs[0]);
  var  = *(int *)mxGetData(prhs[1]);
  int i = 1;
  if (var & PCC_TILED)           /* if to use a tiled version */
    tilesize = *(int *)mxGetData(prhs[++i]); /* get tile size */
  if (var & PCC_THREAD)          /* if to use a threaded version */
    nthreads = *(int *)mxGetData(prhs[++i]); /* get number of threads */
  
  ALLOC_MEM_FOR_RES(             /* allocate memory for result */
          ((size_t)N*(size_t)(N-1))/2)
          
  if (var & PCC_COBL) {
    if (var & PCC_THREAD)        /* cobl. tiling / threads */
      rc = pccx(data, res, N, T, var, 0, nthreads);
    else                         /* cobl. tiling / no threads */
      rc = pccx(data, res, N, T, var, 0);
  } else {
    if (var & PCC_TILED) {
      if (var & PCC_THREAD)      /* std. tiling  / threads */
        rc = pccx(data, res, N, T, var, tilesize, nthreads);
      else                       /* std. tiling  / no threads */
        rc = pccx(data, res, N, T, var, tilesize);
    } else {
      if (var & PCC_THREAD)      /* no tiling    / threads */
        rc = pccx(data, res, N, T, var, nthreads);
      else                       /* no tiling    / no threads */
        rc = pccx(data, res, N, T, var);
    }
  }

  if (rc!=0)
    mexErrMsgTxt("Error.\n");
} /* mexFunction() */
