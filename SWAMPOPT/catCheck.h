/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * catCheck.h
 *
 * Code generation for function 'catCheck'
 *
 */

#ifndef CATCHECK_H
#define CATCHECK_H

/* Include files */
#include "get_MILP_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  void b_sparse_catCheck(const emxArray_int32_T *varargin_1_colidx, int
    varargin_1_m, int varargin_1_n, const emxArray_int32_T *varargin_2_colidx,
    int varargin_2_m, int varargin_2_n, const emxArray_int32_T
    *varargin_3_colidx, int varargin_3_m, int varargin_3_n, const
    emxArray_int32_T *varargin_4_colidx, int varargin_4_m, int varargin_4_n,
    const emxArray_int32_T *varargin_5_colidx, int varargin_5_m, int
    varargin_5_n, const d_sparse varargin_6, const d_sparse varargin_7, const
    d_sparse varargin_8, const d_sparse varargin_9, const d_sparse varargin_10,
    const d_sparse varargin_11, const d_sparse varargin_12, const d_sparse
    varargin_13, int *cnnz, int *cnrows, int *cncols);
  void c_sparse_catCheck(const emxArray_int32_T *varargin_1_colidx, int
    varargin_1_m, int varargin_1_n, const emxArray_real_T *varargin_2, const
    emxArray_int32_T *varargin_3_colidx, int varargin_3_m, int varargin_3_n,
    const emxArray_int32_T *varargin_4_colidx, int varargin_4_m, int
    varargin_4_n, const emxArray_int32_T *varargin_5_colidx, int varargin_5_m,
    int varargin_5_n, const emxArray_real_T *varargin_6, const d_sparse
    varargin_7, const d_sparse varargin_8, int *cnnz, int *cnrows, int *cncols);
  void sparse_catCheck(int varargin_1_m, int varargin_1_n, const emxArray_real_T
                       *varargin_2, int varargin_3_m, int varargin_3_n, const
                       emxArray_real_T *varargin_4, int varargin_5_m, int
                       varargin_5_n, int varargin_6_m, int varargin_6_n, int
                       varargin_7_m, int varargin_7_n, const emxArray_real_T
                       *varargin_8, int varargin_9_m, int varargin_9_n, const
                       d_sparse varargin_10, int *cnnz, int *cnrows, int *cncols);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (catCheck.h) */
