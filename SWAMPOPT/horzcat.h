/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * horzcat.h
 *
 * Code generation for function 'horzcat'
 *
 */

#ifndef HORZCAT_H
#define HORZCAT_H

/* Include files */
#include "get_MILP_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  void b_sparse_horzcat(int varargin_1_m, int varargin_1_n, const
                        emxArray_real_T *varargin_2, int varargin_3_m, int
                        varargin_3_n, const emxArray_real_T *varargin_4, int
                        varargin_5_m, int varargin_5_n, int varargin_6_m, int
                        varargin_6_n, int varargin_7_m, int varargin_7_n, const
                        emxArray_real_T *varargin_8, int varargin_9_m, int
                        varargin_9_n, const d_sparse varargin_10, d_sparse *c);
  void c_sparse_horzcat(const emxArray_real_T *varargin_1, const emxArray_real_T
                        *varargin_2, int varargin_3_m, int varargin_3_n, int
                        varargin_4_m, int varargin_4_n, int varargin_5_m, int
                        varargin_5_n, int varargin_6_m, int varargin_6_n, const
                        emxArray_real_T *varargin_7, int varargin_8_m, int
                        varargin_8_n, d_sparse *c);
  void d_sparse_horzcat(int varargin_1_m, int varargin_1_n, int varargin_2_m,
                        int varargin_2_n, const emxArray_real_T *varargin_3, int
                        varargin_4_m, int varargin_4_n, const emxArray_real_T
                        *varargin_5, d_sparse *c);
  void e_sparse_horzcat(int varargin_1_m, int varargin_1_n, int varargin_2_m,
                        int varargin_2_n, const emxArray_real_T *varargin_3, int
                        varargin_4_m, int varargin_4_n, int varargin_5_m, int
                        varargin_5_n, const emxArray_real_T *varargin_6, int
                        varargin_7_m, int varargin_7_n, d_sparse *c);
  void f_sparse_horzcat(const emxArray_real_T *varargin_1_d, const
                        emxArray_int32_T *varargin_1_colidx, const
                        emxArray_int32_T *varargin_1_rowidx, int varargin_1_m,
                        int varargin_1_n, int varargin_2_m, int varargin_2_n,
                        const emxArray_real_T *varargin_3, d_sparse *c);
  void g_sparse_horzcat(const emxArray_real_T *varargin_1, int varargin_2_m, int
                        varargin_2_n, const emxArray_real_T *varargin_3,
                        d_sparse *c);
  void h_sparse_horzcat(int varargin_1_m, int varargin_1_n, const
                        emxArray_real_T *varargin_2, int varargin_3_m, int
                        varargin_3_n, int varargin_4_m, int varargin_4_n,
                        d_sparse *c);
  void i_sparse_horzcat(int varargin_1_m, int varargin_1_n, const
                        emxArray_real_T *varargin_2, int varargin_3_m, int
                        varargin_3_n, const emxArray_real_T *varargin_4, int
                        varargin_5_m, int varargin_5_n, int varargin_6_m, int
                        varargin_6_n, d_sparse *c);
  void j_sparse_horzcat(int varargin_1_m, int varargin_1_n, const
                        emxArray_real_T *varargin_2, const emxArray_real_T
                        *varargin_3, int varargin_4_m, int varargin_4_n, int
                        varargin_5_m, int varargin_5_n, d_sparse *c);
  void k_sparse_horzcat(int varargin_1_m, int varargin_1_n, const
                        emxArray_real_T *varargin_2, const emxArray_real_T
                        *varargin_3_d, const emxArray_int32_T *varargin_3_colidx,
                        const emxArray_int32_T *varargin_3_rowidx, int
                        varargin_3_m, int varargin_3_n, int varargin_4_m, int
                        varargin_4_n, int varargin_5_m, int varargin_5_n,
                        d_sparse *c);
  void l_sparse_horzcat(int varargin_1_m, int varargin_1_n, const
                        emxArray_real_T *varargin_2_d, const emxArray_int32_T
                        *varargin_2_colidx, const emxArray_int32_T
                        *varargin_2_rowidx, int varargin_2_m, int varargin_2_n,
                        const emxArray_real_T *varargin_3_d, const
                        emxArray_int32_T *varargin_3_colidx, const
                        emxArray_int32_T *varargin_3_rowidx, int varargin_3_m,
                        int varargin_3_n, int varargin_4_m, int varargin_4_n,
                        d_sparse *c);
  void m_sparse_horzcat(const emxArray_real_T *varargin_1, int varargin_2_m, int
                        varargin_2_n, d_sparse *c);
  void n_sparse_horzcat(int varargin_1_m, int varargin_1_n, const
                        emxArray_real_T *varargin_2, const emxArray_real_T
                        *varargin_3, int varargin_4_m, int varargin_4_n,
                        d_sparse *c);
  void o_sparse_horzcat(int K,int varargin_1_m, int varargin_1_n, const
                        emxArray_real_T *varargin_2, int varargin_3_m, int
                        varargin_3_n, d_sparse *c);
  void p_sparse_horzcat(int K,int varargin_1_m, int varargin_1_n, const
                        emxArray_real_T *varargin_2_d, const emxArray_int32_T
                        *varargin_2_colidx, const emxArray_int32_T
                        *varargin_2_rowidx, int varargin_2_m, int varargin_2_n,
                        const emxArray_real_T *varargin_3, int varargin_4_m, int
                        varargin_4_n, d_sparse *c);
  void q_sparse_horzcat(int varargin_1_m, int varargin_1_n, const
                        emxArray_real_T *varargin_2, int varargin_3_m, int
                        varargin_3_n, const emxArray_real_T *varargin_4,
                        d_sparse *c);
  void r_sparse_horzcat(const emxArray_real_T *varargin_1_d, const
                        emxArray_int32_T *varargin_1_colidx, const
                        emxArray_int32_T *varargin_1_rowidx, int varargin_1_m,
                        int varargin_1_n, int varargin_2_m, int varargin_2_n,
                        d_sparse *c);
  void sparse_horzcat(const emxArray_real_T *varargin_1, int varargin_2_m, int
                      varargin_2_n, const emxArray_real_T *varargin_3, int
                      varargin_4_m, int varargin_4_n, int varargin_5_m, int
                      varargin_5_n, int varargin_6_m, int varargin_6_n, const
                      emxArray_real_T *varargin_7, int varargin_8_m, int
                      varargin_8_n, int varargin_9_m, int varargin_9_n, d_sparse
                      *c);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (horzcat.h) */
