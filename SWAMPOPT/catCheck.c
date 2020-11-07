/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * catCheck.c
 *
 * Code generation for function 'catCheck'
 *
 */

/* Include files */
#include "catCheck.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_sparse_catCheck(const emxArray_int32_T *varargin_1_colidx, int
  varargin_1_m, int varargin_1_n, const emxArray_int32_T *varargin_2_colidx, int
  varargin_2_m, int varargin_2_n, const emxArray_int32_T *varargin_3_colidx, int
  varargin_3_m, int varargin_3_n, const emxArray_int32_T *varargin_4_colidx, int
  varargin_4_m, int varargin_4_n, const emxArray_int32_T *varargin_5_colidx, int
  varargin_5_m, int varargin_5_n, const d_sparse varargin_6, const d_sparse
  varargin_7, const d_sparse varargin_8, const d_sparse varargin_9, const
  d_sparse varargin_10, const d_sparse varargin_11, const d_sparse varargin_12,
  const d_sparse varargin_13, int *cnnz, int *cnrows, int *cncols)
{
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T c_isAcceptableEmpty_tmp;
  boolean_T d_isAcceptableEmpty_tmp;
  boolean_T e_isAcceptableEmpty_tmp;
  boolean_T f_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T g_isAcceptableEmpty_tmp;
  boolean_T h_isAcceptableEmpty_tmp;
  boolean_T i_isAcceptableEmpty_tmp;
  boolean_T isAcceptableEmpty_tmp;
  boolean_T j_isAcceptableEmpty_tmp;
  boolean_T k_isAcceptableEmpty_tmp;
  boolean_T l_isAcceptableEmpty_tmp;
  boolean_T m_isAcceptableEmpty_tmp;
  *cncols = varargin_1_n;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  foundSize = !isAcceptableEmpty_tmp;
  if ((varargin_2_m == 0) || (varargin_2_n == 0)) {
    b_isAcceptableEmpty_tmp = true;
  } else {
    b_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      *cncols = varargin_2_n;
    }
  }

  if ((varargin_3_m == 0) || (varargin_3_n == 0)) {
    c_isAcceptableEmpty_tmp = true;
  } else {
    c_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      *cncols = varargin_3_n;
    }
  }

  if ((varargin_4_m == 0) || (varargin_4_n == 0)) {
    d_isAcceptableEmpty_tmp = true;
  } else {
    d_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      *cncols = varargin_4_n;
    }
  }

  if ((varargin_5_m == 0) || (varargin_5_n == 0)) {
    e_isAcceptableEmpty_tmp = true;
  } else {
    e_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      *cncols = varargin_5_n;
    }
  }

  if ((varargin_6.m == 0) || (varargin_6.n == 0)) {
    f_isAcceptableEmpty_tmp = true;
  } else {
    f_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      *cncols = varargin_6.n;
    }
  }

  if ((varargin_7.m == 0) || (varargin_7.n == 0)) {
    g_isAcceptableEmpty_tmp = true;
  } else {
    g_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      *cncols = varargin_7.n;
    }
  }

  if ((varargin_8.m == 0) || (varargin_8.n == 0)) {
    h_isAcceptableEmpty_tmp = true;
  } else {
    h_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      *cncols = varargin_8.n;
    }
  }

  if ((varargin_9.m == 0) || (varargin_9.n == 0)) {
    i_isAcceptableEmpty_tmp = true;
  } else {
    i_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      *cncols = varargin_9.n;
    }
  }

  if ((varargin_10.m == 0) || (varargin_10.n == 0)) {
    j_isAcceptableEmpty_tmp = true;
  } else {
    j_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      *cncols = varargin_10.n;
    }
  }

  if ((varargin_11.m == 0) || (varargin_11.n == 0)) {
    k_isAcceptableEmpty_tmp = true;
  } else {
    k_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      *cncols = varargin_11.n;
    }
  }

  if ((varargin_12.m == 0) || (varargin_12.n == 0)) {
    l_isAcceptableEmpty_tmp = true;
  } else {
    l_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      *cncols = varargin_12.n;
    }
  }

  if ((varargin_13.m == 0) || (varargin_13.n == 0)) {
    m_isAcceptableEmpty_tmp = true;
  } else {
    m_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (isAcceptableEmpty_tmp && b_isAcceptableEmpty_tmp &&
              c_isAcceptableEmpty_tmp && d_isAcceptableEmpty_tmp &&
              e_isAcceptableEmpty_tmp && f_isAcceptableEmpty_tmp &&
              g_isAcceptableEmpty_tmp && h_isAcceptableEmpty_tmp &&
              i_isAcceptableEmpty_tmp && j_isAcceptableEmpty_tmp &&
              k_isAcceptableEmpty_tmp && l_isAcceptableEmpty_tmp &&
              m_isAcceptableEmpty_tmp);
  if ((!m_isAcceptableEmpty_tmp) && (!foundSize)) {
    *cncols = varargin_13.n;
  }

  *cnnz = 0;
  *cnrows = 0;
  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    *cnnz = varargin_1_colidx->data[varargin_1_colidx->size[0] - 1] - 1;
    *cnrows = varargin_1_m;
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_2_colidx->data[varargin_2_colidx->size[0] - 1]) -
      1;
    *cnrows += varargin_2_m;
  }

  if (allEmpty || (!c_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_3_colidx->data[varargin_3_colidx->size[0] - 1]) -
      1;
    *cnrows += varargin_3_m;
  }

  if (allEmpty || (!d_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_4_colidx->data[varargin_4_colidx->size[0] - 1]) -
      1;
    *cnrows += varargin_4_m;
  }

  if (allEmpty || (!e_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_5_colidx->data[varargin_5_colidx->size[0] - 1]) -
      1;
    *cnrows += varargin_5_m;
  }

  if (allEmpty || (!f_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_6.colidx->data[varargin_6.colidx->size[0] - 1]) -
      1;
    *cnrows += varargin_6.m;
  }

  if (allEmpty || (!g_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_7.colidx->data[varargin_7.colidx->size[0] - 1]) -
      1;
    *cnrows += varargin_7.m;
  }

  if (allEmpty || (!h_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_8.colidx->data[varargin_8.colidx->size[0] - 1]) -
      1;
    *cnrows += varargin_8.m;
  }

  if (allEmpty || (!i_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_9.colidx->data[varargin_9.colidx->size[0] - 1]) -
      1;
    *cnrows += varargin_9.m;
  }

  if (allEmpty || (!j_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_10.colidx->data[varargin_10.colidx->size[0] - 1])
      - 1;
    *cnrows += varargin_10.m;
  }

  if (allEmpty || (!k_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_11.colidx->data[varargin_11.colidx->size[0] - 1])
      - 1;
    *cnrows += varargin_11.m;
  }

  if (allEmpty || (!l_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_12.colidx->data[varargin_12.colidx->size[0] - 1])
      - 1;
    *cnrows += varargin_12.m;
  }

  if (allEmpty || (!m_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_13.colidx->data[varargin_13.colidx->size[0] - 1])
      - 1;
    *cnrows += varargin_13.m;
  }
}

void c_sparse_catCheck(const emxArray_int32_T *varargin_1_colidx, int
  varargin_1_m, int varargin_1_n, const emxArray_real_T *varargin_2, const
  emxArray_int32_T *varargin_3_colidx, int varargin_3_m, int varargin_3_n, const
  emxArray_int32_T *varargin_4_colidx, int varargin_4_m, int varargin_4_n, const
  emxArray_int32_T *varargin_5_colidx, int varargin_5_m, int varargin_5_n, const
  emxArray_real_T *varargin_6, const d_sparse varargin_7, const d_sparse
  varargin_8, int *cnnz, int *cnrows, int *cncols)
{
  int i;
  int k;
  int n;
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T c_isAcceptableEmpty_tmp;
  boolean_T d_isAcceptableEmpty_tmp;
  boolean_T e_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  *cncols = varargin_1_n;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  foundSize = !isAcceptableEmpty_tmp;
  isAcceptableEmpty = ((varargin_2->size[0] == 0) || (varargin_2->size[1] == 0));
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    *cncols = varargin_2->size[1];
  }

  if ((varargin_3_m == 0) || (varargin_3_n == 0)) {
    b_isAcceptableEmpty_tmp = true;
  } else {
    b_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      *cncols = varargin_3_n;
    }
  }

  if ((varargin_4_m == 0) || (varargin_4_n == 0)) {
    c_isAcceptableEmpty_tmp = true;
  } else {
    c_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      *cncols = varargin_4_n;
    }
  }

  if ((varargin_5_m == 0) || (varargin_5_n == 0)) {
    d_isAcceptableEmpty_tmp = true;
  } else {
    d_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (isAcceptableEmpty_tmp && isAcceptableEmpty &&
              b_isAcceptableEmpty_tmp && c_isAcceptableEmpty_tmp &&
              d_isAcceptableEmpty_tmp);
  if ((!d_isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    *cncols = varargin_5_n;
  }

  isAcceptableEmpty = ((varargin_6->size[0] == 0) || (varargin_6->size[1] == 0));
  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    *cncols = varargin_6->size[1];
  }

  if ((varargin_7.m == 0) || (varargin_7.n == 0)) {
    isAcceptableEmpty = true;
  } else {
    isAcceptableEmpty = false;
  }

  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    *cncols = varargin_7.n;
  }

  if ((varargin_8.m == 0) || (varargin_8.n == 0)) {
    e_isAcceptableEmpty_tmp = true;
  } else {
    e_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && e_isAcceptableEmpty_tmp);
  if ((!e_isAcceptableEmpty_tmp) && (!foundSize)) {
    *cncols = varargin_8.n;
  }

  *cnnz = 0;
  *cnrows = 0;
  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    *cnnz = varargin_1_colidx->data[varargin_1_colidx->size[0] - 1] - 1;
    *cnrows = varargin_1_m;
  }

  if (allEmpty || ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0))) {
    n = 0;
    i = varargin_2->size[0] * varargin_2->size[1];
    for (k = 0; k < i; k++) {
      if (varargin_2->data[k] != 0.0) {
        n++;
      }
    }

    *cnnz += n;
    *cnrows += varargin_2->size[0];
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_3_colidx->data[varargin_3_colidx->size[0] - 1]) -
      1;
    *cnrows += varargin_3_m;
  }

  if (allEmpty || (!c_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_4_colidx->data[varargin_4_colidx->size[0] - 1]) -
      1;
    *cnrows += varargin_4_m;
  }

  if (allEmpty || (!d_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_5_colidx->data[varargin_5_colidx->size[0] - 1]) -
      1;
    *cnrows += varargin_5_m;
  }

  if (allEmpty || ((varargin_6->size[0] != 0) && (varargin_6->size[1] != 0))) {
    n = 0;
    i = varargin_6->size[0] * varargin_6->size[1];
    for (k = 0; k < i; k++) {
      if (varargin_6->data[k] != 0.0) {
        n++;
      }
    }

    *cnnz += n;
    *cnrows += varargin_6->size[0];
  }

  if (allEmpty || (!isAcceptableEmpty)) {
    *cnnz = (*cnnz + varargin_7.colidx->data[varargin_7.colidx->size[0] - 1]) -
      1;
    *cnrows += varargin_7.m;
  }

  if (allEmpty || (!e_isAcceptableEmpty_tmp)) {
    *cnnz = (*cnnz + varargin_8.colidx->data[varargin_8.colidx->size[0] - 1]) -
      1;
    *cnrows += varargin_8.m;
  }
}

void sparse_catCheck(int varargin_1_m, int varargin_1_n, const emxArray_real_T
                     *varargin_2, int varargin_3_m, int varargin_3_n, const
                     emxArray_real_T *varargin_4, int varargin_5_m, int
                     varargin_5_n, int varargin_6_m, int varargin_6_n, int
                     varargin_7_m, int varargin_7_n, const emxArray_real_T
                     *varargin_8, int varargin_9_m, int varargin_9_n, const
                     d_sparse varargin_10, int *cnnz, int *cnrows, int *cncols)
{
  int i;
  int k;
  int n;
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T c_isAcceptableEmpty_tmp;
  boolean_T d_isAcceptableEmpty_tmp;
  boolean_T e_isAcceptableEmpty_tmp;
  boolean_T f_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  *cnrows = varargin_1_m;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  foundSize = !isAcceptableEmpty_tmp;
  isAcceptableEmpty = ((varargin_2->size[0] == 0) || (varargin_2->size[1] == 0));
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    *cnrows = varargin_2->size[0];
  }

  if ((varargin_3_m == 0) || (varargin_3_n == 0)) {
    b_isAcceptableEmpty_tmp = true;
  } else {
    b_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (isAcceptableEmpty_tmp && isAcceptableEmpty &&
              b_isAcceptableEmpty_tmp);
  if ((!b_isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    *cnrows = varargin_3_m;
  }

  isAcceptableEmpty = ((varargin_4->size[0] == 0) || (varargin_4->size[1] == 0));
  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    *cnrows = varargin_4->size[0];
  }

  if ((varargin_5_m == 0) || (varargin_5_n == 0)) {
    c_isAcceptableEmpty_tmp = true;
  } else {
    c_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && c_isAcceptableEmpty_tmp);
  if ((!c_isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    *cnrows = varargin_5_m;
  }

  if ((varargin_6_m == 0) || (varargin_6_n == 0)) {
    d_isAcceptableEmpty_tmp = true;
  } else {
    d_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && d_isAcceptableEmpty_tmp);
  if ((!d_isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    *cnrows = varargin_6_m;
  }

  if ((varargin_7_m == 0) || (varargin_7_n == 0)) {
    e_isAcceptableEmpty_tmp = true;
  } else {
    e_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && e_isAcceptableEmpty_tmp);
  if ((!e_isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    *cnrows = varargin_7_m;
  }

  isAcceptableEmpty = ((varargin_8->size[0] == 0) || (varargin_8->size[1] == 0));
  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    *cnrows = varargin_8->size[0];
  }

  if ((varargin_9_m == 0) || (varargin_9_n == 0)) {
    isAcceptableEmpty = true;
  } else {
    isAcceptableEmpty = false;
  }

  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    *cnrows = varargin_9_m;
  }

  if ((varargin_10.m == 0) || (varargin_10.n == 0)) {
    f_isAcceptableEmpty_tmp = true;
  } else {
    f_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && f_isAcceptableEmpty_tmp);
  if ((!f_isAcceptableEmpty_tmp) && (!foundSize)) {
    *cnrows = varargin_10.m;
  }

  *cnnz = 0;
  *cncols = 0;
  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    *cncols = varargin_1_n;
  }

  if (allEmpty || ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0))) {
    *cnnz = 0;
    i = varargin_2->size[0] * varargin_2->size[1];
    for (k = 0; k < i; k++) {
      if (varargin_2->data[k] != 0.0) {
        (*cnnz)++;
      }
    }

    *cncols += varargin_2->size[1];
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    *cncols += varargin_3_n;
  }

  if (allEmpty || ((varargin_4->size[0] != 0) && (varargin_4->size[1] != 0))) {
    n = 0;
    i = varargin_4->size[0] * varargin_4->size[1];
    for (k = 0; k < i; k++) {
      if (varargin_4->data[k] != 0.0) {
        n++;
      }
    }

    *cnnz += n;
    *cncols += varargin_4->size[1];
  }

  if (allEmpty || (!c_isAcceptableEmpty_tmp)) {
    *cncols += varargin_5_n;
  }

  if (allEmpty || (!d_isAcceptableEmpty_tmp)) {
    *cncols += varargin_6_n;
  }

  if (allEmpty || (!e_isAcceptableEmpty_tmp)) {
    *cncols += varargin_7_n;
  }

  if (allEmpty || ((varargin_8->size[0] != 0) && (varargin_8->size[1] != 0))) {
    n = 0;
    i = varargin_8->size[0] * varargin_8->size[1];
    for (k = 0; k < i; k++) {
      if (varargin_8->data[k] != 0.0) {
        n++;
      }
    }

    *cnnz += n;
    *cncols += varargin_8->size[1];
  }

  if (allEmpty || (!isAcceptableEmpty)) {
    *cncols += varargin_9_n;
  }

  if (allEmpty || (!f_isAcceptableEmpty_tmp)) {
    *cncols += varargin_10.n;
  }
}

/* End of code generation (catCheck.c) */
