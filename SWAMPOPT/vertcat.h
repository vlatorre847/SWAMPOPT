/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * vertcat.h
 *
 * Code generation for function 'vertcat'
 *
 */

#ifndef VERTCAT_H
#define VERTCAT_H

/* Include files */
#include "get_MILP_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  void b_sparse_vertcat(const emxArray_real_T *varargin_1_d, const
                        emxArray_int32_T *varargin_1_colidx, const
                        emxArray_int32_T *varargin_1_rowidx, int varargin_1_m,
                        int varargin_1_n, const emxArray_real_T *varargin_2_d,
                        const emxArray_int32_T *varargin_2_colidx, const
                        emxArray_int32_T *varargin_2_rowidx, int varargin_2_m,
                        int varargin_2_n, const emxArray_real_T *varargin_3_d,
                        const emxArray_int32_T *varargin_3_colidx, const
                        emxArray_int32_T *varargin_3_rowidx, int varargin_3_m,
                        int varargin_3_n, d_sparse *c);
  void c_sparse_vertcat(const emxArray_real_T *varargin_1_d, const
                        emxArray_int32_T *varargin_1_colidx, const
                        emxArray_int32_T *varargin_1_rowidx, int varargin_1_m,
                        int varargin_1_n, const emxArray_real_T *varargin_2_d,
                        const emxArray_int32_T *varargin_2_colidx, const
                        emxArray_int32_T *varargin_2_rowidx, int varargin_2_m,
                        int varargin_2_n, const emxArray_real_T *varargin_3_d,
                        const emxArray_int32_T *varargin_3_colidx, const
                        emxArray_int32_T *varargin_3_rowidx, int varargin_3_m,
                        int varargin_3_n, const d_sparse varargin_4, const
                        d_sparse varargin_5, const d_sparse varargin_6, const
                        d_sparse varargin_7, const d_sparse varargin_8, const
                        d_sparse varargin_9, const d_sparse varargin_10, const
                        d_sparse varargin_11, const d_sparse varargin_12, const
                        d_sparse varargin_13, d_sparse *c);
  void d_sparse_vertcat(const emxArray_real_T *varargin_1_d, const
                        emxArray_int32_T *varargin_1_colidx, const
                        emxArray_int32_T *varargin_1_rowidx, int varargin_1_m,
                        int varargin_1_n, const emxArray_real_T *varargin_2,
                        const emxArray_real_T *varargin_3_d, const
                        emxArray_int32_T *varargin_3_colidx, const
                        emxArray_int32_T *varargin_3_rowidx, int varargin_3_m,
                        int varargin_3_n, const emxArray_real_T *varargin_4_d,
                        const emxArray_int32_T *varargin_4_colidx, const
                        emxArray_int32_T *varargin_4_rowidx, int varargin_4_m,
                        int varargin_4_n, const d_sparse varargin_5, const
                        emxArray_real_T *varargin_6, const d_sparse varargin_7,
                        const d_sparse varargin_8, d_sparse *c);
  void sparse_vertcat(const emxArray_real_T *varargin_1_d, const
                      emxArray_int32_T *varargin_1_colidx, const
                      emxArray_int32_T *varargin_1_rowidx, int varargin_1_m, int
                      varargin_1_n, const emxArray_real_T *varargin_2_d, const
                      emxArray_int32_T *varargin_2_colidx, const
                      emxArray_int32_T *varargin_2_rowidx, int varargin_2_m, int
                      varargin_2_n, d_sparse *c);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (vertcat.h) */
