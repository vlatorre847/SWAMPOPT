/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * vertcat.c
 *
 * Code generation for function 'vertcat'
 *
 */

/* Include files */
#include "vertcat.h"
#include "catCheck.h"
#include "fillIn.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_sparse_vertcat(const emxArray_real_T *varargin_1_d, const
                      emxArray_int32_T *varargin_1_colidx, const
                      emxArray_int32_T *varargin_1_rowidx, int varargin_1_m, int
                      varargin_1_n, const emxArray_real_T *varargin_2_d, const
                      emxArray_int32_T *varargin_2_colidx, const
                      emxArray_int32_T *varargin_2_rowidx, int varargin_2_m, int
                      varargin_2_n, const emxArray_real_T *varargin_3_d, const
                      emxArray_int32_T *varargin_3_colidx, const
                      emxArray_int32_T *varargin_3_rowidx, int varargin_3_m, int
                      varargin_3_n, d_sparse *c)
{
  int ccol;
  int cnfixeddim;
  int cnvardim;
  int i;
  int kp;
  int kpend;
  int kpend_tmp;
  int kpstart;
  int numalloc;
  boolean_T allEmpty;
  boolean_T emptyflag_idx_0;
  boolean_T emptyflag_idx_2;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty_tmp;
  cnfixeddim = varargin_1_n;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    emptyflag_idx_0 = true;
  } else {
    emptyflag_idx_0 = false;
  }

  foundSize = !emptyflag_idx_0;
  if ((varargin_2_m == 0) || (varargin_2_n == 0)) {
    emptyflag_idx_2 = true;
  } else {
    emptyflag_idx_2 = false;
    if (!foundSize) {
      foundSize = true;
      cnfixeddim = varargin_2_n;
    }
  }

  if ((varargin_3_m == 0) || (varargin_3_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  allEmpty = (emptyflag_idx_0 && emptyflag_idx_2 && isAcceptableEmpty_tmp);
  if ((!isAcceptableEmpty_tmp) && (!foundSize)) {
    cnfixeddim = varargin_3_n;
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || (!emptyflag_idx_0)) {
    numalloc = varargin_1_colidx->data[varargin_1_colidx->size[0] - 1] - 1;
    cnvardim = varargin_1_m;
  }

  if (allEmpty || (!emptyflag_idx_2)) {
    numalloc = (numalloc + varargin_2_colidx->data[varargin_2_colidx->size[0] -
                1]) - 1;
    cnvardim += varargin_2_m;
  }

  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    numalloc = (numalloc + varargin_3_colidx->data[varargin_3_colidx->size[0] -
                1]) - 1;
    cnvardim += varargin_3_m;
  }

  c->m = cnvardim;
  c->n = cnfixeddim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  i = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, i);
  for (i = 0; i < numalloc; i++) {
    c->d->data[i] = 0.0;
  }

  c->maxnz = numalloc;
  i = c->colidx->size[0];
  c->colidx->size[0] = cnfixeddim + 1;
  emxEnsureCapacity_int32_T(c->colidx, i);
  c->colidx->data[0] = 1;
  i = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, i);
  for (i = 0; i < numalloc; i++) {
    c->rowidx->data[i] = 0;
  }

  for (cnvardim = 0; cnvardim < cnfixeddim; cnvardim++) {
    c->colidx->data[cnvardim + 1] = 1;
  }

  sparse_fillIn(c);
  cnvardim = -1;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    emptyflag_idx_0 = true;
  } else {
    emptyflag_idx_0 = false;
  }

  if ((varargin_2_m == 0) || (varargin_2_n == 0)) {
    foundSize = true;
  } else {
    foundSize = false;
  }

  if ((varargin_3_m == 0) || (varargin_3_n == 0)) {
    emptyflag_idx_2 = true;
  } else {
    emptyflag_idx_2 = false;
  }

  i = c->n;
  for (ccol = 0; ccol < i; ccol++) {
    numalloc = 0;
    if (!emptyflag_idx_0) {
      cnfixeddim = cnvardim;
      kpstart = varargin_1_colidx->data[ccol];
      kpend_tmp = varargin_1_colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cnfixeddim++;
        c->rowidx->data[cnfixeddim] = varargin_1_rowidx->data[kp - 1];
        c->d->data[cnfixeddim] = varargin_1_d->data[kp - 1];
      }

      cnvardim = (cnvardim + kpend_tmp) - varargin_1_colidx->data[ccol];
      numalloc = varargin_1_m;
    }

    if (!foundSize) {
      cnfixeddim = cnvardim;
      kpstart = varargin_2_colidx->data[ccol];
      kpend_tmp = varargin_2_colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cnfixeddim++;
        c->rowidx->data[cnfixeddim] = varargin_2_rowidx->data[kp - 1] + numalloc;
        c->d->data[cnfixeddim] = varargin_2_d->data[kp - 1];
      }

      cnvardim = (cnvardim + kpend_tmp) - varargin_2_colidx->data[ccol];
      numalloc += varargin_2_m;
    }

    if (!emptyflag_idx_2) {
      cnfixeddim = cnvardim;
      kpstart = varargin_3_colidx->data[ccol];
      kpend_tmp = varargin_3_colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cnfixeddim++;
        c->rowidx->data[cnfixeddim] = varargin_3_rowidx->data[kp - 1] + numalloc;
        c->d->data[cnfixeddim] = varargin_3_d->data[kp - 1];
      }

      cnvardim = (cnvardim + kpend_tmp) - varargin_3_colidx->data[ccol];
    }

    c->colidx->data[ccol + 1] = cnvardim + 2;
  }
}

void c_sparse_vertcat(const emxArray_real_T *varargin_1_d, const
                      emxArray_int32_T *varargin_1_colidx, const
                      emxArray_int32_T *varargin_1_rowidx, int varargin_1_m, int
                      varargin_1_n, const emxArray_real_T *varargin_2_d, const
                      emxArray_int32_T *varargin_2_colidx, const
                      emxArray_int32_T *varargin_2_rowidx, int varargin_2_m, int
                      varargin_2_n, const emxArray_real_T *varargin_3_d, const
                      emxArray_int32_T *varargin_3_colidx, const
                      emxArray_int32_T *varargin_3_rowidx, int varargin_3_m, int
                      varargin_3_n, const d_sparse varargin_4, const d_sparse
                      varargin_5, const d_sparse varargin_6, const d_sparse
                      varargin_7, const d_sparse varargin_8, const d_sparse
                      varargin_9, const d_sparse varargin_10, const d_sparse
                      varargin_11, const d_sparse varargin_12, const d_sparse
                      varargin_13, d_sparse *c)
{
  int ccol;
  int cidx;
  int crowoffs;
  int i;
  int kp;
  int kpend;
  int kpend_tmp;
  int kpstart;
  int numalloc;
  boolean_T emptyflag_idx_0;
  boolean_T emptyflag_idx_1;
  boolean_T emptyflag_idx_10;
  boolean_T emptyflag_idx_11;
  boolean_T emptyflag_idx_12;
  boolean_T emptyflag_idx_2;
  boolean_T emptyflag_idx_3;
  boolean_T emptyflag_idx_4;
  boolean_T emptyflag_idx_5;
  boolean_T emptyflag_idx_6;
  boolean_T emptyflag_idx_7;
  boolean_T emptyflag_idx_8;
  boolean_T emptyflag_idx_9;
  b_sparse_catCheck(varargin_1_colidx, varargin_1_m, varargin_1_n,
                    varargin_2_colidx, varargin_2_m, varargin_2_n,
                    varargin_3_colidx, varargin_3_m, varargin_3_n,
                    varargin_4.colidx, varargin_4.m, varargin_4.n,
                    varargin_5.colidx, varargin_5.m, varargin_5.n, varargin_6,
                    varargin_7, varargin_8, varargin_9, varargin_10, varargin_11,
                    varargin_12, varargin_13, &numalloc, &c->m, &crowoffs);
  c->n = crowoffs;
  if (numalloc < 1) {
    numalloc = 1;
  }

  i = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, i);
  for (i = 0; i < numalloc; i++) {
    c->d->data[i] = 0.0;
  }

  c->maxnz = numalloc;
  i = c->colidx->size[0];
  c->colidx->size[0] = crowoffs + 1;
  emxEnsureCapacity_int32_T(c->colidx, i);
  c->colidx->data[0] = 1;
  i = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, i);
  for (i = 0; i < numalloc; i++) {
    c->rowidx->data[i] = 0;
  }

  for (numalloc = 0; numalloc < crowoffs; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  numalloc = -1;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    emptyflag_idx_0 = true;
  } else {
    emptyflag_idx_0 = false;
  }

  if ((varargin_2_m == 0) || (varargin_2_n == 0)) {
    emptyflag_idx_1 = true;
  } else {
    emptyflag_idx_1 = false;
  }

  if ((varargin_3_m == 0) || (varargin_3_n == 0)) {
    emptyflag_idx_2 = true;
  } else {
    emptyflag_idx_2 = false;
  }

  if ((varargin_4.m == 0) || (varargin_4.n == 0)) {
    emptyflag_idx_3 = true;
  } else {
    emptyflag_idx_3 = false;
  }

  if ((varargin_5.m == 0) || (varargin_5.n == 0)) {
    emptyflag_idx_4 = true;
  } else {
    emptyflag_idx_4 = false;
  }

  if ((varargin_6.m == 0) || (varargin_6.n == 0)) {
    emptyflag_idx_5 = true;
  } else {
    emptyflag_idx_5 = false;
  }

  if ((varargin_7.m == 0) || (varargin_7.n == 0)) {
    emptyflag_idx_6 = true;
  } else {
    emptyflag_idx_6 = false;
  }

  if ((varargin_8.m == 0) || (varargin_8.n == 0)) {
    emptyflag_idx_7 = true;
  } else {
    emptyflag_idx_7 = false;
  }

  if ((varargin_9.m == 0) || (varargin_9.n == 0)) {
    emptyflag_idx_8 = true;
  } else {
    emptyflag_idx_8 = false;
  }

  if ((varargin_10.m == 0) || (varargin_10.n == 0)) {
    emptyflag_idx_9 = true;
  } else {
    emptyflag_idx_9 = false;
  }

  if ((varargin_11.m == 0) || (varargin_11.n == 0)) {
    emptyflag_idx_10 = true;
  } else {
    emptyflag_idx_10 = false;
  }

  if ((varargin_12.m == 0) || (varargin_12.n == 0)) {
    emptyflag_idx_11 = true;
  } else {
    emptyflag_idx_11 = false;
  }

  if ((varargin_13.m == 0) || (varargin_13.n == 0)) {
    emptyflag_idx_12 = true;
  } else {
    emptyflag_idx_12 = false;
  }

  i = c->n;
  for (ccol = 0; ccol < i; ccol++) {
    crowoffs = 0;
    if (!emptyflag_idx_0) {
      cidx = numalloc;
      kpstart = varargin_1_colidx->data[ccol];
      kpend_tmp = varargin_1_colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_1_rowidx->data[kp - 1];
        c->d->data[cidx] = varargin_1_d->data[kp - 1];
      }

      numalloc = (numalloc + kpend_tmp) - varargin_1_colidx->data[ccol];
      crowoffs = varargin_1_m;
    }

    if (!emptyflag_idx_1) {
      cidx = numalloc;
      kpstart = varargin_2_colidx->data[ccol];
      kpend_tmp = varargin_2_colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_2_rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_2_d->data[kp - 1];
      }

      numalloc = (numalloc + kpend_tmp) - varargin_2_colidx->data[ccol];
      crowoffs += varargin_2_m;
    }

    if (!emptyflag_idx_2) {
      cidx = numalloc;
      kpstart = varargin_3_colidx->data[ccol];
      kpend_tmp = varargin_3_colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_3_rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_3_d->data[kp - 1];
      }

      numalloc = (numalloc + kpend_tmp) - varargin_3_colidx->data[ccol];
      crowoffs += varargin_3_m;
    }

    if (!emptyflag_idx_3) {
      cidx = numalloc;
      kpstart = varargin_4.colidx->data[ccol];
      kpend_tmp = varargin_4.colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_4.rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_4.d->data[kp - 1];
      }

      numalloc = (numalloc + kpend_tmp) - varargin_4.colidx->data[ccol];
      crowoffs += varargin_4.m;
    }

    if (!emptyflag_idx_4) {
      cidx = numalloc;
      kpstart = varargin_5.colidx->data[ccol];
      kpend_tmp = varargin_5.colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_5.rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_5.d->data[kp - 1];
      }

      numalloc = (numalloc + kpend_tmp) - varargin_5.colidx->data[ccol];
      crowoffs += varargin_5.m;
    }

    if (!emptyflag_idx_5) {
      cidx = numalloc;
      kpstart = varargin_6.colidx->data[ccol];
      kpend_tmp = varargin_6.colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_6.rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_6.d->data[kp - 1];
      }

      numalloc = (numalloc + kpend_tmp) - varargin_6.colidx->data[ccol];
      crowoffs += varargin_6.m;
    }

    if (!emptyflag_idx_6) {
      cidx = numalloc;
      kpstart = varargin_7.colidx->data[ccol];
      kpend_tmp = varargin_7.colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_7.rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_7.d->data[kp - 1];
      }

      numalloc = (numalloc + kpend_tmp) - varargin_7.colidx->data[ccol];
      crowoffs += varargin_7.m;
    }

    if (!emptyflag_idx_7) {
      cidx = numalloc;
      kpstart = varargin_8.colidx->data[ccol];
      kpend_tmp = varargin_8.colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_8.rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_8.d->data[kp - 1];
      }

      numalloc = (numalloc + kpend_tmp) - varargin_8.colidx->data[ccol];
      crowoffs += varargin_8.m;
    }

    if (!emptyflag_idx_8) {
      cidx = numalloc;
      kpstart = varargin_9.colidx->data[ccol];
      kpend_tmp = varargin_9.colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_9.rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_9.d->data[kp - 1];
      }

      numalloc = (numalloc + kpend_tmp) - varargin_9.colidx->data[ccol];
      crowoffs += varargin_9.m;
    }

    if (!emptyflag_idx_9) {
      cidx = numalloc;
      kpstart = varargin_10.colidx->data[ccol];
      kpend_tmp = varargin_10.colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_10.rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_10.d->data[kp - 1];
      }

      numalloc = (numalloc + kpend_tmp) - varargin_10.colidx->data[ccol];
      crowoffs += varargin_10.m;
    }

    if (!emptyflag_idx_10) {
      cidx = numalloc;
      kpstart = varargin_11.colidx->data[ccol];
      kpend_tmp = varargin_11.colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_11.rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_11.d->data[kp - 1];
      }

      numalloc = (numalloc + kpend_tmp) - varargin_11.colidx->data[ccol];
      crowoffs += varargin_11.m;
    }

    if (!emptyflag_idx_11) {
      cidx = numalloc;
      kpstart = varargin_12.colidx->data[ccol];
      kpend_tmp = varargin_12.colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_12.rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_12.d->data[kp - 1];
      }

      numalloc = (numalloc + kpend_tmp) - varargin_12.colidx->data[ccol];
      crowoffs += varargin_12.m;
    }

    if (!emptyflag_idx_12) {
      cidx = numalloc;
      kpstart = varargin_13.colidx->data[ccol];
      kpend_tmp = varargin_13.colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_13.rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_13.d->data[kp - 1];
      }

      numalloc = (numalloc + kpend_tmp) - varargin_13.colidx->data[ccol];
    }

    c->colidx->data[ccol + 1] = numalloc + 2;
  }
}

void d_sparse_vertcat(const emxArray_real_T *varargin_1_d, const
                      emxArray_int32_T *varargin_1_colidx, const
                      emxArray_int32_T *varargin_1_rowidx, int varargin_1_m, int
                      varargin_1_n, const emxArray_real_T *varargin_2, const
                      emxArray_real_T *varargin_3_d, const emxArray_int32_T
                      *varargin_3_colidx, const emxArray_int32_T
                      *varargin_3_rowidx, int varargin_3_m, int varargin_3_n,
                      const emxArray_real_T *varargin_4_d, const
                      emxArray_int32_T *varargin_4_colidx, const
                      emxArray_int32_T *varargin_4_rowidx, int varargin_4_m, int
                      varargin_4_n, const d_sparse varargin_5, const
                      emxArray_real_T *varargin_6, const d_sparse varargin_7,
                      const d_sparse varargin_8, d_sparse *c)
{
  int ccol;
  int cidx;
  int cncols;
  int crowoffs;
  int i;
  int kp;
  int kpend;
  int numalloc;
  int nzCount;
  boolean_T emptyflag_idx_0;
  boolean_T emptyflag_idx_1;
  boolean_T emptyflag_idx_2;
  boolean_T emptyflag_idx_3;
  boolean_T emptyflag_idx_4;
  boolean_T emptyflag_idx_5;
  boolean_T emptyflag_idx_6;
  boolean_T emptyflag_idx_7;
  c_sparse_catCheck(varargin_1_colidx, varargin_1_m, varargin_1_n, varargin_2,
                    varargin_3_colidx, varargin_3_m, varargin_3_n,
                    varargin_4_colidx, varargin_4_m, varargin_4_n,
                    varargin_5.colidx, varargin_5.m, varargin_5.n, varargin_6,
                    varargin_7, varargin_8, &numalloc, &c->m, &cncols);
  c->n = cncols;
  if (numalloc < 1) {
    numalloc = 1;
  }

  i = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, i);
  for (i = 0; i < numalloc; i++) {
    c->d->data[i] = 0.0;
  }

  c->maxnz = numalloc;
  i = c->colidx->size[0];
  c->colidx->size[0] = cncols + 1;
  emxEnsureCapacity_int32_T(c->colidx, i);
  c->colidx->data[0] = 1;
  i = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, i);
  for (i = 0; i < numalloc; i++) {
    c->rowidx->data[i] = 0;
  }

  for (numalloc = 0; numalloc < cncols; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  nzCount = -1;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    emptyflag_idx_0 = true;
  } else {
    emptyflag_idx_0 = false;
  }

  emptyflag_idx_1 = ((varargin_2->size[0] == 0) || (varargin_2->size[1] == 0));
  if ((varargin_3_m == 0) || (varargin_3_n == 0)) {
    emptyflag_idx_2 = true;
  } else {
    emptyflag_idx_2 = false;
  }

  if ((varargin_4_m == 0) || (varargin_4_n == 0)) {
    emptyflag_idx_3 = true;
  } else {
    emptyflag_idx_3 = false;
  }

  if ((varargin_5.m == 0) || (varargin_5.n == 0)) {
    emptyflag_idx_4 = true;
  } else {
    emptyflag_idx_4 = false;
  }

  emptyflag_idx_5 = ((varargin_6->size[0] == 0) || (varargin_6->size[1] == 0));
  if ((varargin_7.m == 0) || (varargin_7.n == 0)) {
    emptyflag_idx_6 = true;
  } else {
    emptyflag_idx_6 = false;
  }

  if ((varargin_8.m == 0) || (varargin_8.n == 0)) {
    emptyflag_idx_7 = true;
  } else {
    emptyflag_idx_7 = false;
  }

  i = c->n;
  for (ccol = 0; ccol < i; ccol++) {
    crowoffs = 0;
    if (!emptyflag_idx_0) {
      cidx = nzCount;
      numalloc = varargin_1_colidx->data[ccol];
      cncols = varargin_1_colidx->data[ccol + 1];
      kpend = cncols - 1;
      for (kp = numalloc; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_1_rowidx->data[kp - 1];
        c->d->data[cidx] = varargin_1_d->data[kp - 1];
      }

      nzCount = (nzCount + cncols) - varargin_1_colidx->data[ccol];
      crowoffs = varargin_1_m;
    }

    if (!emptyflag_idx_1) {
      numalloc = varargin_2->size[0];
      cidx = nzCount;
      for (cncols = 0; cncols < numalloc; cncols++) {
        if (varargin_2->data[cncols + varargin_2->size[0] * ccol] != 0.0) {
          cidx++;
          c->rowidx->data[cidx] = (cncols + crowoffs) + 1;
          c->d->data[cidx] = 1.0;
        }
      }

      nzCount = cidx;
      crowoffs += varargin_2->size[0];
    }

    if (!emptyflag_idx_2) {
      cidx = nzCount;
      numalloc = varargin_3_colidx->data[ccol];
      cncols = varargin_3_colidx->data[ccol + 1];
      kpend = cncols - 1;
      for (kp = numalloc; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_3_rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_3_d->data[kp - 1];
      }

      nzCount = (nzCount + cncols) - varargin_3_colidx->data[ccol];
      crowoffs += varargin_3_m;
    }

    if (!emptyflag_idx_3) {
      cidx = nzCount;
      numalloc = varargin_4_colidx->data[ccol];
      cncols = varargin_4_colidx->data[ccol + 1];
      kpend = cncols - 1;
      for (kp = numalloc; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_4_rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_4_d->data[kp - 1];
      }

      nzCount = (nzCount + cncols) - varargin_4_colidx->data[ccol];
      crowoffs += varargin_4_m;
    }

    if (!emptyflag_idx_4) {
      cidx = nzCount;
      numalloc = varargin_5.colidx->data[ccol];
      cncols = varargin_5.colidx->data[ccol + 1];
      kpend = cncols - 1;
      for (kp = numalloc; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_5.rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_5.d->data[kp - 1];
      }

      nzCount = (nzCount + cncols) - varargin_5.colidx->data[ccol];
      crowoffs += varargin_5.m;
    }

    if (!emptyflag_idx_5) {
      numalloc = varargin_6->size[0];
      cidx = nzCount;
      for (cncols = 0; cncols < numalloc; cncols++) {
        if (varargin_6->data[cncols + varargin_6->size[0] * ccol] != 0.0) {
          cidx++;
          c->rowidx->data[cidx] = (cncols + crowoffs) + 1;
          c->d->data[cidx] = 1.0;
        }
      }

      nzCount = cidx;
      crowoffs += varargin_6->size[0];
    }

    if (!emptyflag_idx_6) {
      cidx = nzCount;
      numalloc = varargin_7.colidx->data[ccol];
      cncols = varargin_7.colidx->data[ccol + 1];
      kpend = cncols - 1;
      for (kp = numalloc; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_7.rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_7.d->data[kp - 1];
      }

      nzCount = (nzCount + cncols) - varargin_7.colidx->data[ccol];
      crowoffs += varargin_7.m;
    }

    if (!emptyflag_idx_7) {
      cidx = nzCount;
      numalloc = varargin_8.colidx->data[ccol];
      cncols = varargin_8.colidx->data[ccol + 1];
      kpend = cncols - 1;
      for (kp = numalloc; kp <= kpend; kp++) {
        cidx++;
        c->rowidx->data[cidx] = varargin_8.rowidx->data[kp - 1] + crowoffs;
        c->d->data[cidx] = varargin_8.d->data[kp - 1];
      }

      nzCount = (nzCount + cncols) - varargin_8.colidx->data[ccol];
    }

    c->colidx->data[ccol + 1] = nzCount + 2;
  }
}

void sparse_vertcat(const emxArray_real_T *varargin_1_d, const emxArray_int32_T *
                    varargin_1_colidx, const emxArray_int32_T *varargin_1_rowidx,
                    int varargin_1_m, int varargin_1_n, const emxArray_real_T
                    *varargin_2_d, const emxArray_int32_T *varargin_2_colidx,
                    const emxArray_int32_T *varargin_2_rowidx, int varargin_2_m,
                    int varargin_2_n, d_sparse *c)
{
  int ccol;
  int cnfixeddim;
  int cnvardim;
  int i;
  int kp;
  int kpend;
  int kpend_tmp;
  int kpstart;
  int numalloc;
  boolean_T allEmpty;
  boolean_T emptyflag_idx_0;
  boolean_T emptyflag_idx_1;
  cnfixeddim = varargin_1_n;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    emptyflag_idx_0 = true;
  } else {
    emptyflag_idx_0 = false;
  }

  if ((varargin_2_m == 0) || (varargin_2_n == 0)) {
    emptyflag_idx_1 = true;
  } else {
    emptyflag_idx_1 = false;
  }

  allEmpty = (emptyflag_idx_0 && emptyflag_idx_1);
  if ((!emptyflag_idx_1) && emptyflag_idx_0) {
    cnfixeddim = varargin_2_n;
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || (!emptyflag_idx_0)) {
    numalloc = varargin_1_colidx->data[varargin_1_colidx->size[0] - 1] - 1;
    cnvardim = varargin_1_m;
  }

  if (allEmpty || (!emptyflag_idx_1)) {
    numalloc = (numalloc + varargin_2_colidx->data[varargin_2_colidx->size[0] -
                1]) - 1;
    cnvardim += varargin_2_m;
  }

  c->m = cnvardim;
  c->n = cnfixeddim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  i = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, i);
  for (i = 0; i < numalloc; i++) {
    c->d->data[i] = 0.0;
  }

  c->maxnz = numalloc;
  i = c->colidx->size[0];
  c->colidx->size[0] = cnfixeddim + 1;
  emxEnsureCapacity_int32_T(c->colidx, i);
  c->colidx->data[0] = 1;
  i = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, i);
  for (i = 0; i < numalloc; i++) {
    c->rowidx->data[i] = 0;
  }

  for (cnvardim = 0; cnvardim < cnfixeddim; cnvardim++) {
    c->colidx->data[cnvardim + 1] = 1;
  }

  sparse_fillIn(c);
  cnvardim = -1;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    emptyflag_idx_0 = true;
  } else {
    emptyflag_idx_0 = false;
  }

  if ((varargin_2_m == 0) || (varargin_2_n == 0)) {
    emptyflag_idx_1 = true;
  } else {
    emptyflag_idx_1 = false;
  }

  i = c->n;
  for (ccol = 0; ccol < i; ccol++) {
    numalloc = 0;
    if (!emptyflag_idx_0) {
      cnfixeddim = cnvardim;
      kpstart = varargin_1_colidx->data[ccol];
      kpend_tmp = varargin_1_colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cnfixeddim++;
        c->rowidx->data[cnfixeddim] = varargin_1_rowidx->data[kp - 1];
        c->d->data[cnfixeddim] = varargin_1_d->data[kp - 1];
      }

      cnvardim = (cnvardim + kpend_tmp) - varargin_1_colidx->data[ccol];
      numalloc = varargin_1_m;
    }

    if (!emptyflag_idx_1) {
      cnfixeddim = cnvardim;
      kpstart = varargin_2_colidx->data[ccol];
      kpend_tmp = varargin_2_colidx->data[ccol + 1];
      kpend = kpend_tmp - 1;
      for (kp = kpstart; kp <= kpend; kp++) {
        cnfixeddim++;
        c->rowidx->data[cnfixeddim] = varargin_2_rowidx->data[kp - 1] + numalloc;
        c->d->data[cnfixeddim] = varargin_2_d->data[kp - 1];
      }

      cnvardim = (cnvardim + kpend_tmp) - varargin_2_colidx->data[ccol];
    }

    c->colidx->data[ccol + 1] = cnvardim + 2;
  }
}

/* End of code generation (vertcat.c) */
