/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sparse1.h
 *
 * Code generation for function 'sparse1'
 *
 */

#ifndef SPARSE1_H
#define SPARSE1_H

/* Include files */
#include "get_MILP_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  void sparse_ctranspose(const emxArray_real_T *this_d, const emxArray_int32_T
    *this_colidx, const emxArray_int32_T *this_rowidx, int this_n,
    emxArray_real_T *y_d, emxArray_int32_T *y_colidx, emxArray_int32_T *y_rowidx,
    int *y_m);
  void sparse_full(const emxArray_real_T *this_d, const emxArray_int32_T
                   *this_colidx, const emxArray_int32_T *this_rowidx, int this_m,
                   emxArray_real_T *y);
  void sparse_ne(const emxArray_real_T *a_d, const emxArray_int32_T *a_colidx,
                 const emxArray_int32_T *a_rowidx, int a_m, int a_n,
                 emxArray_boolean_T *s_d, emxArray_int32_T *s_colidx,
                 emxArray_int32_T *s_rowidx, int *s_m, int *s_n);
  void sparse_nonzeros(const emxArray_real_T *this_d, const emxArray_int32_T
                       *this_colidx, emxArray_real_T *y);
  void sparse_parenAssign(d_sparse *this, int rhs_n, double varargin_1);
  void sparse_parenReference(const emxArray_real_T *this_d, const
    emxArray_int32_T *this_colidx, const emxArray_int32_T *this_rowidx, int
    this_n, const emxArray_real_T *varargin_1, emxArray_real_T *s_d,
    emxArray_int32_T *s_colidx, emxArray_int32_T *s_rowidx, int *s_m, int *s_n,
    int *s_maxnz);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (sparse1.h) */
