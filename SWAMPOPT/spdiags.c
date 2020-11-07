/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * spdiags.c
 *
 * Code generation for function 'spdiags'
 *
 */

/* Include files */
#include "spdiags.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"
#include "sparse.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void b_spdiags(int K,const double arg1[K], double arg3, double arg4, emxArray_real_T *
               res1_d, emxArray_int32_T *res1_colidx, emxArray_int32_T
               *res1_rowidx, int *res1_m, int *res1_n)
{
  d_sparse expl_temp;
  emxArray_int32_T *aCols;
  emxArray_int32_T *aRows;
  emxArray_int32_T *aRows_tmp;
  emxArray_int32_T *r;
  emxArray_real_T *aDat;
  emxArray_real_T *idx;
  emxArray_real_T *y;
  double len[2];
  double minAdjustedDim_data_idx_0;
  int b_loop_ub;
  int c_loop_ub;
  int i;
  int k;
  int loop_ub;
  int trueCount;
  trueCount = 0;
  if ((0.0 >= -arg3 + 1.0) && (0.0 <= arg4 - 1.0)) {
    trueCount = 1;
  }

  len[0] = 0.0;
  if (0 <= trueCount - 1) {
    if ((arg3 < arg4) || rtIsNaN(arg4)) {
      minAdjustedDim_data_idx_0 = arg3;
    } else {
      minAdjustedDim_data_idx_0 = arg4;
    }

    len[1] = (minAdjustedDim_data_idx_0 - 1.0) + 1.0;
  }

  emxInit_int32_T(&aRows, 1);
  loop_ub = (int)len[trueCount];
  i = aRows->size[0];
  aRows->size[0] = loop_ub;
  emxEnsureCapacity_int32_T(aRows, i);
  for (i = 0; i < loop_ub; i++) {
    aRows->data[i] = 0;
  }

  emxInit_int32_T(&aCols, 1);
  i = aCols->size[0];
  aCols->size[0] = loop_ub;
  emxEnsureCapacity_int32_T(aCols, i);
  for (i = 0; i < loop_ub; i++) {
    aCols->data[i] = 0;
  }

  emxInit_real_T(&aDat, 1);
  i = aDat->size[0];
  aDat->size[0] = loop_ub;
  emxEnsureCapacity_real_T(aDat, i);
  for (i = 0; i < loop_ub; i++) {
    aDat->data[i] = 0.0;
  }

  emxInit_real_T(&idx, 2);
  emxInit_real_T(&y, 2);
  if (0 <= trueCount - 1) {
    if (rtIsNaN(minAdjustedDim_data_idx_0)) {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity_real_T(y, i);
      y->data[0] = rtNaN;
    } else if (minAdjustedDim_data_idx_0 < 1.0) {
      y->size[0] = 1;
      y->size[1] = 0;
    } else if (rtIsInf(minAdjustedDim_data_idx_0) && (1.0 ==
                minAdjustedDim_data_idx_0)) {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity_real_T(y, i);
      y->data[0] = rtNaN;
    } else {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      loop_ub = (int)floor(minAdjustedDim_data_idx_0 - 1.0);
      y->size[1] = loop_ub + 1;
      emxEnsureCapacity_real_T(y, i);
      for (i = 0; i <= loop_ub; i++) {
        y->data[i] = (double)i + 1.0;
      }
    }

    if (rtIsNaN(len[1])) {
      i = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      idx->size[1] = 1;
      emxEnsureCapacity_real_T(idx, i);
      idx->data[0] = rtNaN;
    } else if (len[1] < 1.0) {
      idx->size[0] = 1;
      idx->size[1] = 0;
    } else if (rtIsInf(len[1]) && (1.0 == len[1])) {
      i = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      idx->size[1] = 1;
      emxEnsureCapacity_real_T(idx, i);
      idx->data[0] = rtNaN;
    } else {
      i = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      loop_ub = (int)floor(len[1] - 1.0);
      idx->size[1] = loop_ub + 1;
      emxEnsureCapacity_real_T(idx, i);
      for (i = 0; i <= loop_ub; i++) {
        idx->data[i] = (double)i + 1.0;
      }
    }

    b_loop_ub = idx->size[0] * idx->size[1];
    c_loop_ub = y->size[1];
  }

  emxInit_int32_T(&r, 2);
  emxInit_int32_T(&aRows_tmp, 1);
  for (k = 0; k < trueCount; k++) {
    i = r->size[0] * r->size[1];
    r->size[0] = 1;
    r->size[1] = idx->size[1];
    emxEnsureCapacity_int32_T(r, i);
    for (i = 0; i < b_loop_ub; i++) {
      r->data[i] = (int)idx->data[i];
    }

    i = aRows_tmp->size[0];
    aRows_tmp->size[0] = y->size[1];
    emxEnsureCapacity_int32_T(aRows_tmp, i);
    for (i = 0; i < c_loop_ub; i++) {
      aRows_tmp->data[i] = (int)y->data[i];
    }

    loop_ub = r->size[0] * r->size[1];
    for (i = 0; i < loop_ub; i++) {
      aRows->data[r->data[i] - 1] = aRows_tmp->data[i];
    }

    loop_ub = r->size[0] * r->size[1];
    for (i = 0; i < loop_ub; i++) {
      aCols->data[r->data[i] - 1] = aRows_tmp->data[i];
    }

    loop_ub = r->size[1];
    for (i = 0; i < loop_ub; i++) {
      aDat->data[r->data[i] - 1] = arg1[aRows_tmp->data[i] - 1];
    }
  }

  emxFree_real_T(&y);
  emxFree_int32_T(&aRows_tmp);
  emxFree_int32_T(&r);
  emxFree_real_T(&idx);
  emxInitStruct_sparse(&expl_temp);
  b_sparse(aRows, aCols, aDat, arg3, arg4, &expl_temp);
  i = res1_d->size[0];
  res1_d->size[0] = expl_temp.d->size[0];
  emxEnsureCapacity_real_T(res1_d, i);
  loop_ub = expl_temp.d->size[0];
  emxFree_real_T(&aDat);
  emxFree_int32_T(&aCols);
  emxFree_int32_T(&aRows);
  for (i = 0; i < loop_ub; i++) {
    res1_d->data[i] = expl_temp.d->data[i];
  }

  i = res1_colidx->size[0];
  res1_colidx->size[0] = expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(res1_colidx, i);
  loop_ub = expl_temp.colidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    res1_colidx->data[i] = expl_temp.colidx->data[i];
  }

  i = res1_rowidx->size[0];
  res1_rowidx->size[0] = expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(res1_rowidx, i);
  loop_ub = expl_temp.rowidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    res1_rowidx->data[i] = expl_temp.rowidx->data[i];
  }

  *res1_m = expl_temp.m;
  *res1_n = expl_temp.n;
  emxFreeStruct_sparse(&expl_temp);
}

void spdiags(const emxArray_real_T *arg1, double arg3, double arg4,
             emxArray_real_T *res1_d, emxArray_int32_T *res1_colidx,
             emxArray_int32_T *res1_rowidx, int *res1_m, int *res1_n)
{
  d_sparse expl_temp;
  emxArray_int32_T *aCols;
  emxArray_int32_T *aRows;
  emxArray_int32_T *aRows_tmp;
  emxArray_int32_T *r;
  emxArray_real_T *aDat;
  emxArray_real_T *idx;
  emxArray_real_T *y;
  double len[2];
  double minAdjustedDim_data_idx_0;
  int b_loop_ub;
  int c_loop_ub;
  int i;
  int k;
  int loop_ub;
  int trueCount;
  trueCount = 0;
  if ((0.0 >= -arg3 + 1.0) && (0.0 <= arg4 - 1.0)) {
    trueCount = 1;
  }

  len[0] = 0.0;
  if (0 <= trueCount - 1) {
    if ((arg3 < arg4) || rtIsNaN(arg4)) {
      minAdjustedDim_data_idx_0 = arg3;
    } else {
      minAdjustedDim_data_idx_0 = arg4;
    }

    len[1] = (minAdjustedDim_data_idx_0 - 1.0) + 1.0;
  }

  emxInit_int32_T(&aRows, 1);
  loop_ub = (int)len[trueCount];
  i = aRows->size[0];
  aRows->size[0] = loop_ub;
  emxEnsureCapacity_int32_T(aRows, i);
  for (i = 0; i < loop_ub; i++) {
    aRows->data[i] = 0;
  }

  emxInit_int32_T(&aCols, 1);
  i = aCols->size[0];
  aCols->size[0] = loop_ub;
  emxEnsureCapacity_int32_T(aCols, i);
  for (i = 0; i < loop_ub; i++) {
    aCols->data[i] = 0;
  }

  emxInit_real_T(&aDat, 1);
  i = aDat->size[0];
  aDat->size[0] = loop_ub;
  emxEnsureCapacity_real_T(aDat, i);
  for (i = 0; i < loop_ub; i++) {
    aDat->data[i] = 0.0;
  }

  emxInit_real_T(&idx, 2);
  emxInit_real_T(&y, 2);
  if (0 <= trueCount - 1) {
    if (rtIsNaN(minAdjustedDim_data_idx_0)) {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity_real_T(y, i);
      y->data[0] = rtNaN;
    } else if (minAdjustedDim_data_idx_0 < 1.0) {
      y->size[0] = 1;
      y->size[1] = 0;
    } else if (rtIsInf(minAdjustedDim_data_idx_0) && (1.0 ==
                minAdjustedDim_data_idx_0)) {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity_real_T(y, i);
      y->data[0] = rtNaN;
    } else {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      loop_ub = (int)floor(minAdjustedDim_data_idx_0 - 1.0);
      y->size[1] = loop_ub + 1;
      emxEnsureCapacity_real_T(y, i);
      for (i = 0; i <= loop_ub; i++) {
        y->data[i] = (double)i + 1.0;
      }
    }

    if (rtIsNaN(len[1])) {
      i = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      idx->size[1] = 1;
      emxEnsureCapacity_real_T(idx, i);
      idx->data[0] = rtNaN;
    } else if (len[1] < 1.0) {
      idx->size[0] = 1;
      idx->size[1] = 0;
    } else if (rtIsInf(len[1]) && (1.0 == len[1])) {
      i = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      idx->size[1] = 1;
      emxEnsureCapacity_real_T(idx, i);
      idx->data[0] = rtNaN;
    } else {
      i = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      loop_ub = (int)floor(len[1] - 1.0);
      idx->size[1] = loop_ub + 1;
      emxEnsureCapacity_real_T(idx, i);
      for (i = 0; i <= loop_ub; i++) {
        idx->data[i] = (double)i + 1.0;
      }
    }

    b_loop_ub = idx->size[0] * idx->size[1];
    c_loop_ub = y->size[1];
  }

  emxInit_int32_T(&r, 2);
  emxInit_int32_T(&aRows_tmp, 1);
  for (k = 0; k < trueCount; k++) {
    i = r->size[0] * r->size[1];
    r->size[0] = 1;
    r->size[1] = idx->size[1];
    emxEnsureCapacity_int32_T(r, i);
    for (i = 0; i < b_loop_ub; i++) {
      r->data[i] = (int)idx->data[i];
    }

    i = aRows_tmp->size[0];
    aRows_tmp->size[0] = y->size[1];
    emxEnsureCapacity_int32_T(aRows_tmp, i);
    for (i = 0; i < c_loop_ub; i++) {
      aRows_tmp->data[i] = (int)y->data[i];
    }

    loop_ub = r->size[0] * r->size[1];
    for (i = 0; i < loop_ub; i++) {
      aRows->data[r->data[i] - 1] = aRows_tmp->data[i];
    }

    loop_ub = r->size[0] * r->size[1];
    for (i = 0; i < loop_ub; i++) {
      aCols->data[r->data[i] - 1] = aRows_tmp->data[i];
    }

    loop_ub = r->size[0] * r->size[1];
    for (i = 0; i < loop_ub; i++) {
      aDat->data[r->data[i] - 1] = arg1->data[aRows_tmp->data[i] - 1];
    }
  }

  emxFree_real_T(&y);
  emxFree_int32_T(&aRows_tmp);
  emxFree_int32_T(&r);
  emxFree_real_T(&idx);
  emxInitStruct_sparse(&expl_temp);
  b_sparse(aRows, aCols, aDat, arg3, arg4, &expl_temp);
  i = res1_d->size[0];
  res1_d->size[0] = expl_temp.d->size[0];
  emxEnsureCapacity_real_T(res1_d, i);
  loop_ub = expl_temp.d->size[0];
  emxFree_real_T(&aDat);
  emxFree_int32_T(&aCols);
  emxFree_int32_T(&aRows);
  for (i = 0; i < loop_ub; i++) {
    res1_d->data[i] = expl_temp.d->data[i];
  }

  i = res1_colidx->size[0];
  res1_colidx->size[0] = expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(res1_colidx, i);
  loop_ub = expl_temp.colidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    res1_colidx->data[i] = expl_temp.colidx->data[i];
  }

  i = res1_rowidx->size[0];
  res1_rowidx->size[0] = expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(res1_rowidx, i);
  loop_ub = expl_temp.rowidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    res1_rowidx->data[i] = expl_temp.rowidx->data[i];
  }

  *res1_m = expl_temp.m;
  *res1_n = expl_temp.n;
  emxFreeStruct_sparse(&expl_temp);
}

/* End of code generation (spdiags.c) */
