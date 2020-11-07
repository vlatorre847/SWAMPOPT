/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sparse.h
 *
 * Code generation for function 'sparse'
 *
 */

#ifndef SPARSE_H
#define SPARSE_H

/* Include files */
#include "get_MILP_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  void b_sparse(const emxArray_int32_T *varargin_1, const emxArray_int32_T
                *varargin_2, const emxArray_real_T *varargin_3, double
                varargin_4, double varargin_5, d_sparse *y);
  void c_sparse(double varargin_2, emxArray_real_T *y_d, emxArray_int32_T
                *y_colidx, emxArray_int32_T *y_rowidx, int *y_n);
  void sparse(double varargin_1, double varargin_2, emxArray_real_T *y_d,
              emxArray_int32_T *y_colidx, emxArray_int32_T *y_rowidx, int *y_m,
              int *y_n, int *y_maxnz);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (sparse.h) */
