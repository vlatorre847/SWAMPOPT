/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * repmat1.h
 *
 * Code generation for function 'repmat1'
 *
 */

#ifndef REPMAT1_H
#define REPMAT1_H

/* Include files */
#include "get_MILP_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  void sparse_repmat(const emxArray_real_T *A_d, const emxArray_int32_T
                     *A_colidx, const emxArray_int32_T *A_rowidx, int A_m, int
                     A_n, double varargin_2, d_sparse *B);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (repmat1.h) */
