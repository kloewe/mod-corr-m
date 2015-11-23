/*----------------------------------------------------------------------------
  File    : mxCorrAlloc.h
  Author  : Kristian Loewe
----------------------------------------------------------------------------*/
#ifndef MX_CORRALLOC_H
#define MX_CORRALLOC_H

#ifndef REAL
#define REAL float
#endif

#define float  1    /* to check definitions */
#define double 2

#if REAL==float
#define ALLOC_MEM_FOR_RES(s) { plhs[0] = mxCreateNumericMatrix( s, 1, \
      mxSINGLE_CLASS, mxREAL); res = (float *)mxGetData(plhs[0]); }
#elif REAL==double
#define ALLOC_MEM_FOR_RES(s) { plhs[0] = mxCreateNumericMatrix(s, 1, \
      mxDOUBLE_CLASS, mxREAL); res = mxGetPr(plhs[0]); }
#endif

#undef float
#undef double

#endif  /* #ifndef MX_CORRALLOC_H */
