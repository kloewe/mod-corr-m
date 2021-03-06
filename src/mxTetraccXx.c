/*----------------------------------------------------------------------------
  File    : mxTetraccXx.c
  Contents: MEX wrapper for the C-based functions defined in tetracc.h
  Author  : Kristian Loewe
----------------------------------------------------------------------------*/
#include "mex.h"
#include "tetracc.h"
#include "mxCorrAlloc.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{                                 /* --- gateway routine */
  int  N,                         /* number of nodes */
       T,                         /* number of points in time */
       var;                       /* implementation variant */
  int  tile = 0;                  /* tile size */
  int  nthd = 0;                  /* number of threads */
  REAL *data;                     /* input data (nodal time series) */
  REAL *res;                      /* result */
  int  rc = 1;                    /* return condition */

  N = (int)mxGetN(prhs[0]);       /* retrieve input arguments from MATLAB */
  T = (int)mxGetM(prhs[0]);
  data = (REAL *)mxGetData(prhs[0]);
  var  = *(int *)mxGetData(prhs[1]);
  int i = 1;
  if (var & TCC_TILED)            /* if to use tiling */
    tile = *(int *)mxGetData(prhs[++i]); /* get tile size */
  if (var & TCC_THREAD)           /* if to use multi-threading */
    nthd = *(int *)mxGetData(prhs[++i]); /* get number of threads */

  ALLOC_MEM_FOR_RES(              /* allocate memory for result */
          ((size_t)N*(size_t)(N-1))/2)

  if (var & TCC_TILED)            /* if to use tiling */
    if (var & TCC_THREAD)         /*   if to use multi-threaded version */
      rc = tetraccx(data, res, N, T, var, tile, nthd); 
    else                          /*   else: use single-threaded version */
      rc = tetraccx(data, res, N, T, var, tile);
  else                            /* else: no tiling */
    if (var & TCC_THREAD)         /*   if to use multi-threaded version */
      rc = tetraccx(data, res, N, T, var, nthd);
    else                          /*   else: use single-threaded version */
      rc = tetraccx(data, res, N, T, var);

  if (rc)
    mexErrMsgTxt("An error occured.\n");
}  /* mexFunction() */
