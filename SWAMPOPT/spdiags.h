/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * spdiags.h
 *
 * Code generation for function 'spdiags'
 *
 */

#ifndef SPDIAGS_H
#define SPDIAGS_H

/* Include files */
#include "get_MILP_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  void b_spdiags(int K,const double arg1[K], double arg3, double arg4,
                 emxArray_real_T *res1_d, emxArray_int32_T *res1_colidx,
                 emxArray_int32_T *res1_rowidx, int *res1_m, int *res1_n);
  void spdiags(const emxArray_real_T *arg1, double arg3, double arg4,
               emxArray_real_T *res1_d, emxArray_int32_T *res1_colidx,
               emxArray_int32_T *res1_rowidx, int *res1_m, int *res1_n);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (spdiags.h) */
