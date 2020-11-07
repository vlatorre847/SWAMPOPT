/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * horzcat.c
 *
 * Code generation for function 'horzcat'
 *
 */

/* Include files */
#include "horzcat.h"
#include "catCheck.h"
#include "fillIn.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_sparse_horzcat(int varargin_1_m, int varargin_1_n, const emxArray_real_T *
                      varargin_2, int varargin_3_m, int varargin_3_n, const
                      emxArray_real_T *varargin_4, int varargin_5_m, int
                      varargin_5_n, int varargin_6_m, int varargin_6_n, int
                      varargin_7_m, int varargin_7_n, const emxArray_real_T
                      *varargin_8, int varargin_9_m, int varargin_9_n, const
                      d_sparse varargin_10, d_sparse *c)
{
  double dk;
  int ccolidx;
  int cidx;
  int col;
  int ncolk;
  int numalloc;
  int row;
  sparse_catCheck(varargin_1_m, varargin_1_n, varargin_2, varargin_3_m,
                  varargin_3_n, varargin_4, varargin_5_m, varargin_5_n,
                  varargin_6_m, varargin_6_n, varargin_7_m, varargin_7_n,
                  varargin_8, varargin_9_m, varargin_9_n, varargin_10, &numalloc,
                  &c->m, &ncolk);
  c->n = ncolk;
  if (numalloc < 1) {
    numalloc = 1;
  }

  row = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, row);
  for (row = 0; row < numalloc; row++) {
    c->d->data[row] = 0.0;
  }

  c->maxnz = numalloc;
  row = c->colidx->size[0];
  c->colidx->size[0] = ncolk + 1;
  emxEnsureCapacity_int32_T(c->colidx, row);
  c->colidx->data[0] = 1;
  row = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, row);
  for (row = 0; row < numalloc; row++) {
    c->rowidx->data[row] = 0;
  }

  for (numalloc = 0; numalloc < ncolk; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cidx = 1;
  ccolidx = 0;
  if ((varargin_1_m != 0) && (varargin_1_n != 0)) {
    for (col = 0; col < varargin_1_n; col++) {
      ccolidx++;
      c->colidx->data[ccolidx] = 1;
    }
  }

  if ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0)) {
    numalloc = varargin_2->size[0];
    ncolk = varargin_2->size[1];
    cidx = 1;
    for (col = 0; col < ncolk; col++) {
      for (row = 0; row < numalloc; row++) {
        if (varargin_2->data[row + varargin_2->size[0] * col] != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = -1.0;
          cidx++;
        }
      }

      ccolidx++;
      c->colidx->data[ccolidx] = cidx;
    }
  }

  if ((varargin_3_m != 0) && (varargin_3_n != 0)) {
    for (col = 0; col < varargin_3_n; col++) {
      ccolidx++;
      c->colidx->data[ccolidx] = cidx;
    }
  }

  if ((varargin_4->size[0] != 0) && (varargin_4->size[1] != 0)) {
    numalloc = varargin_4->size[0];
    ncolk = varargin_4->size[1];
    for (col = 0; col < ncolk; col++) {
      for (row = 0; row < numalloc; row++) {
        dk = varargin_4->data[row + varargin_4->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = dk;
          cidx++;
        }
      }

      ccolidx++;
      c->colidx->data[ccolidx] = cidx;
    }
  }

  if ((varargin_5_m != 0) && (varargin_5_n != 0)) {
    for (col = 0; col < varargin_5_n; col++) {
      ccolidx++;
      c->colidx->data[ccolidx] = cidx;
    }
  }

  if ((varargin_6_m != 0) && (varargin_6_n != 0)) {
    for (col = 0; col < varargin_6_n; col++) {
      ccolidx++;
      c->colidx->data[ccolidx] = cidx;
    }
  }

  if ((varargin_7_m != 0) && (varargin_7_n != 0)) {
    for (col = 0; col < varargin_7_n; col++) {
      ccolidx++;
      c->colidx->data[ccolidx] = cidx;
    }
  }

  if ((varargin_8->size[0] != 0) && (varargin_8->size[1] != 0)) {
    numalloc = varargin_8->size[0];
    ncolk = varargin_8->size[1];
    for (col = 0; col < ncolk; col++) {
      for (row = 0; row < numalloc; row++) {
        dk = varargin_8->data[row + varargin_8->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = dk;
          cidx++;
        }
      }

      ccolidx++;
      c->colidx->data[ccolidx] = cidx;
    }
  }

  if ((varargin_9_m != 0) && (varargin_9_n != 0)) {
    for (col = 0; col < varargin_9_n; col++) {
      ccolidx++;
      c->colidx->data[ccolidx] = cidx;
    }
  }

  if ((varargin_10.m != 0) && (varargin_10.n != 0)) {
    row = varargin_10.n;
    for (col = 0; col < row; col++) {
      ccolidx++;
      c->colidx->data[ccolidx] = cidx;
    }
  }
}

void c_sparse_horzcat(const emxArray_real_T *varargin_1, const emxArray_real_T
                      *varargin_2, int varargin_3_m, int varargin_3_n, int
                      varargin_4_m, int varargin_4_n, int varargin_5_m, int
                      varargin_5_n, int varargin_6_m, int varargin_6_n, const
                      emxArray_real_T *varargin_7, int varargin_8_m, int
                      varargin_8_n, d_sparse *c)
{
  double dk;
  int cidx;
  int cnvardim;
  int col;
  int nrowk;
  int numalloc;
  int row;
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T c_isAcceptableEmpty_tmp;
  boolean_T d_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1->size[0];
  isAcceptableEmpty = ((varargin_1->size[0] == 0) || (varargin_1->size[1] == 0));
  allEmpty = isAcceptableEmpty;
  foundSize = !isAcceptableEmpty;
  isAcceptableEmpty = ((varargin_2->size[0] == 0) || (varargin_2->size[1] == 0));
  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_2->size[0];
  }

  if ((varargin_3_m == 0) || (varargin_3_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && isAcceptableEmpty_tmp);
  if ((!isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_3_m;
  }

  if ((varargin_4_m == 0) || (varargin_4_n == 0)) {
    b_isAcceptableEmpty_tmp = true;
  } else {
    b_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && b_isAcceptableEmpty_tmp);
  if ((!b_isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_4_m;
  }

  if ((varargin_5_m == 0) || (varargin_5_n == 0)) {
    c_isAcceptableEmpty_tmp = true;
  } else {
    c_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && c_isAcceptableEmpty_tmp);
  if ((!c_isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_5_m;
  }

  if ((varargin_6_m == 0) || (varargin_6_n == 0)) {
    d_isAcceptableEmpty_tmp = true;
  } else {
    d_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && d_isAcceptableEmpty_tmp);
  if ((!d_isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_6_m;
  }

  isAcceptableEmpty = ((varargin_7->size[0] == 0) || (varargin_7->size[1] == 0));
  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_7->size[0];
  }

  if ((varargin_8_m == 0) || (varargin_8_n == 0)) {
    isAcceptableEmpty = true;
  } else {
    isAcceptableEmpty = false;
  }

  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    c->m = varargin_8_m;
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0))) {
    numalloc = 0;
    col = varargin_1->size[0] * varargin_1->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_1->data[cidx] != 0.0) {
        numalloc++;
      }
    }

    cnvardim = varargin_1->size[1];
  }

  if (allEmpty || ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0))) {
    nrowk = 0;
    col = varargin_2->size[0] * varargin_2->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_2->data[cidx] != 0.0) {
        nrowk++;
      }
    }

    numalloc += nrowk;
    cnvardim += varargin_2->size[1];
  }

  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim += varargin_3_n;
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_4_n;
  }

  if (allEmpty || (!c_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_5_n;
  }

  if (allEmpty || (!d_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_6_n;
  }

  if (allEmpty || ((varargin_7->size[0] != 0) && (varargin_7->size[1] != 0))) {
    nrowk = 0;
    col = varargin_7->size[0] * varargin_7->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_7->data[cidx] != 0.0) {
        nrowk++;
      }
    }

    numalloc += nrowk;
    cnvardim += varargin_7->size[1];
  }

  if (allEmpty || (!isAcceptableEmpty)) {
    cnvardim += varargin_8_n;
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  col = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, col);
  for (col = 0; col < numalloc; col++) {
    c->d->data[col] = 0.0;
  }

  c->maxnz = numalloc;
  col = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, col);
  c->colidx->data[0] = 1;
  col = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, col);
  for (col = 0; col < numalloc; col++) {
    c->rowidx->data[col] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cidx = 1;
  numalloc = 0;
  if ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0)) {
    nrowk = varargin_1->size[0];
    cnvardim = varargin_1->size[1];
    cidx = 1;
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        dk = varargin_1->data[row + varargin_1->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0)) {
    nrowk = varargin_2->size[0];
    cnvardim = varargin_2->size[1];
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        dk = varargin_2->data[row + varargin_2->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_3_m != 0) && (varargin_3_n != 0)) {
    for (col = 0; col < varargin_3_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_4_m != 0) && (varargin_4_n != 0)) {
    for (col = 0; col < varargin_4_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_5_m != 0) && (varargin_5_n != 0)) {
    for (col = 0; col < varargin_5_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_6_m != 0) && (varargin_6_n != 0)) {
    for (col = 0; col < varargin_6_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_7->size[0] != 0) && (varargin_7->size[1] != 0)) {
    nrowk = varargin_7->size[0];
    cnvardim = varargin_7->size[1];
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        dk = varargin_7->data[row + varargin_7->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_8_m != 0) && (varargin_8_n != 0)) {
    for (col = 0; col < varargin_8_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }
}

void d_sparse_horzcat(int varargin_1_m, int varargin_1_n, int varargin_2_m, int
                      varargin_2_n, const emxArray_real_T *varargin_3, int
                      varargin_4_m, int varargin_4_n, const emxArray_real_T
                      *varargin_5, d_sparse *c)
{
  double dk;
  int cidx;
  int cnvardim;
  int col;
  int nrowk;
  int numalloc;
  int row;
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T c_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1_m;
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
      c->m = varargin_2_m;
    }
  }

  isAcceptableEmpty = ((varargin_3->size[0] == 0) || (varargin_3->size[1] == 0));
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_3->size[0];
  }

  if ((varargin_4_m == 0) || (varargin_4_n == 0)) {
    c_isAcceptableEmpty_tmp = true;
  } else {
    c_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (isAcceptableEmpty_tmp && b_isAcceptableEmpty_tmp &&
              isAcceptableEmpty && c_isAcceptableEmpty_tmp);
  if ((!c_isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_4_m;
  }

  isAcceptableEmpty = ((varargin_5->size[0] == 0) || (varargin_5->size[1] == 0));
  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    c->m = varargin_5->size[0];
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim = varargin_1_n;
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_2_n;
  }

  if (allEmpty || ((varargin_3->size[0] != 0) && (varargin_3->size[1] != 0))) {
    numalloc = 0;
    col = varargin_3->size[0] * varargin_3->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_3->data[cidx] != 0.0) {
        numalloc++;
      }
    }

    cnvardim += varargin_3->size[1];
  }

  if (allEmpty || (!c_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_4_n;
  }

  if (allEmpty || ((varargin_5->size[0] != 0) && (varargin_5->size[1] != 0))) {
    nrowk = 0;
    col = varargin_5->size[0] * varargin_5->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_5->data[cidx] != 0.0) {
        nrowk++;
      }
    }

    numalloc += nrowk;
    cnvardim += varargin_5->size[1];
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  col = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, col);
  for (col = 0; col < numalloc; col++) {
    c->d->data[col] = 0.0;
  }

  c->maxnz = numalloc;
  col = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, col);
  c->colidx->data[0] = 1;
  col = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, col);
  for (col = 0; col < numalloc; col++) {
    c->rowidx->data[col] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cidx = 0;
  numalloc = 0;
  if ((varargin_1_m != 0) && (varargin_1_n != 0)) {
    for (col = 0; col < varargin_1_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = 1;
    }
  }

  if ((varargin_2_m != 0) && (varargin_2_n != 0)) {
    for (col = 0; col < varargin_2_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = 1;
    }
  }

  if ((varargin_3->size[0] != 0) && (varargin_3->size[1] != 0)) {
    nrowk = varargin_3->size[0];
    cnvardim = varargin_3->size[1];
    cidx = 0;
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        if (varargin_3->data[row + varargin_3->size[0] * col] != 0.0) {
          c->rowidx->data[cidx] = row + 1;
          c->d->data[cidx] = -1.0;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx + 1;
    }
  }

  if ((varargin_4_m != 0) && (varargin_4_n != 0)) {
    for (col = 0; col < varargin_4_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx + 1;
    }
  }

  if ((varargin_5->size[0] != 0) && (varargin_5->size[1] != 0)) {
    nrowk = varargin_5->size[0];
    cnvardim = varargin_5->size[1];
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        dk = varargin_5->data[row + varargin_5->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx] = row + 1;
          c->d->data[cidx] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx + 1;
    }
  }
}

void e_sparse_horzcat(int varargin_1_m, int varargin_1_n, int varargin_2_m, int
                      varargin_2_n, const emxArray_real_T *varargin_3, int
                      varargin_4_m, int varargin_4_n, int varargin_5_m, int
                      varargin_5_n, const emxArray_real_T *varargin_6, int
                      varargin_7_m, int varargin_7_n, d_sparse *c)
{
  double dk;
  int cidx;
  int cnvardim;
  int col;
  int nrowk;
  int numalloc;
  int row;
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T c_isAcceptableEmpty_tmp;
  boolean_T d_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1_m;
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
      c->m = varargin_2_m;
    }
  }

  isAcceptableEmpty = ((varargin_3->size[0] == 0) || (varargin_3->size[1] == 0));
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_3->size[0];
  }

  if ((varargin_4_m == 0) || (varargin_4_n == 0)) {
    c_isAcceptableEmpty_tmp = true;
  } else {
    c_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      c->m = varargin_4_m;
    }
  }

  if ((varargin_5_m == 0) || (varargin_5_n == 0)) {
    d_isAcceptableEmpty_tmp = true;
  } else {
    d_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (isAcceptableEmpty_tmp && b_isAcceptableEmpty_tmp &&
              isAcceptableEmpty && c_isAcceptableEmpty_tmp &&
              d_isAcceptableEmpty_tmp);
  if ((!d_isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_5_m;
  }

  isAcceptableEmpty = ((varargin_6->size[0] == 0) || (varargin_6->size[1] == 0));
  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_6->size[0];
  }

  if ((varargin_7_m == 0) || (varargin_7_n == 0)) {
    isAcceptableEmpty = true;
  } else {
    isAcceptableEmpty = false;
  }

  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    c->m = varargin_7_m;
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim = varargin_1_n;
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_2_n;
  }

  if (allEmpty || ((varargin_3->size[0] != 0) && (varargin_3->size[1] != 0))) {
    numalloc = 0;
    col = varargin_3->size[0] * varargin_3->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_3->data[cidx] != 0.0) {
        numalloc++;
      }
    }

    cnvardim += varargin_3->size[1];
  }

  if (allEmpty || (!c_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_4_n;
  }

  if (allEmpty || (!d_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_5_n;
  }

  if (allEmpty || ((varargin_6->size[0] != 0) && (varargin_6->size[1] != 0))) {
    nrowk = 0;
    col = varargin_6->size[0] * varargin_6->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_6->data[cidx] != 0.0) {
        nrowk++;
      }
    }

    numalloc += nrowk;
    cnvardim += varargin_6->size[1];
  }

  if (allEmpty || (!isAcceptableEmpty)) {
    cnvardim += varargin_7_n;
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  col = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, col);
  for (col = 0; col < numalloc; col++) {
    c->d->data[col] = 0.0;
  }

  c->maxnz = numalloc;
  col = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, col);
  c->colidx->data[0] = 1;
  col = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, col);
  for (col = 0; col < numalloc; col++) {
    c->rowidx->data[col] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cidx = 1;
  numalloc = 0;
  if ((varargin_1_m != 0) && (varargin_1_n != 0)) {
    for (col = 0; col < varargin_1_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = 1;
    }
  }

  if ((varargin_2_m != 0) && (varargin_2_n != 0)) {
    for (col = 0; col < varargin_2_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = 1;
    }
  }

  if ((varargin_3->size[0] != 0) && (varargin_3->size[1] != 0)) {
    nrowk = varargin_3->size[0];
    cnvardim = varargin_3->size[1];
    cidx = 1;
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        if (varargin_3->data[row + varargin_3->size[0] * col] != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = -1.0;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_4_m != 0) && (varargin_4_n != 0)) {
    for (col = 0; col < varargin_4_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_5_m != 0) && (varargin_5_n != 0)) {
    for (col = 0; col < varargin_5_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_6->size[0] != 0) && (varargin_6->size[1] != 0)) {
    nrowk = varargin_6->size[0];
    cnvardim = varargin_6->size[1];
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        dk = varargin_6->data[row + varargin_6->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_7_m != 0) && (varargin_7_n != 0)) {
    for (col = 0; col < varargin_7_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }
}

void f_sparse_horzcat(const emxArray_real_T *varargin_1_d, const
                      emxArray_int32_T *varargin_1_colidx, const
                      emxArray_int32_T *varargin_1_rowidx, int varargin_1_m, int
                      varargin_1_n, int varargin_2_m, int varargin_2_n, const
                      emxArray_real_T *varargin_3, d_sparse *c)
{
  double dk;
  int ccolidx;
  int cnvardim;
  int col;
  int ncolk;
  int numalloc;
  int row;
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1_m;
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
      c->m = varargin_2_m;
    }
  }

  isAcceptableEmpty = ((varargin_3->size[0] == 0) || (varargin_3->size[1] == 0));
  allEmpty = (isAcceptableEmpty_tmp && b_isAcceptableEmpty_tmp &&
              isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    c->m = varargin_3->size[0];
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    numalloc = varargin_1_colidx->data[varargin_1_colidx->size[0] - 1] - 1;
    cnvardim = varargin_1_n;
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_2_n;
  }

  if (allEmpty || ((varargin_3->size[0] != 0) && (varargin_3->size[1] != 0))) {
    ncolk = 0;
    col = varargin_3->size[0] * varargin_3->size[1];
    for (ccolidx = 0; ccolidx < col; ccolidx++) {
      if (varargin_3->data[ccolidx] != 0.0) {
        ncolk++;
      }
    }

    numalloc += ncolk;
    cnvardim += varargin_3->size[1];
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  col = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, col);
  for (col = 0; col < numalloc; col++) {
    c->d->data[col] = 0.0;
  }

  c->maxnz = numalloc;
  col = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, col);
  c->colidx->data[0] = 1;
  col = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, col);
  for (col = 0; col < numalloc; col++) {
    c->rowidx->data[col] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cnvardim = 1;
  ccolidx = 0;
  if ((varargin_1_m != 0) && (varargin_1_n != 0)) {
    cnvardim = -1;
    numalloc = varargin_1_colidx->data[varargin_1_colidx->size[0] - 1];
    for (ncolk = 0; ncolk <= numalloc - 2; ncolk++) {
      cnvardim++;
      c->rowidx->data[cnvardim] = varargin_1_rowidx->data[ncolk];
      c->d->data[cnvardim] = varargin_1_d->data[ncolk];
    }

    for (col = 0; col < varargin_1_n; col++) {
      ccolidx++;
      c->colidx->data[ccolidx] = varargin_1_colidx->data[col + 1];
    }

    cnvardim = varargin_1_colidx->data[varargin_1_colidx->size[0] - 1];
  }

  if ((varargin_2_m != 0) && (varargin_2_n != 0)) {
    for (col = 0; col < varargin_2_n; col++) {
      ccolidx++;
      c->colidx->data[ccolidx] = cnvardim;
    }
  }

  if ((varargin_3->size[0] != 0) && (varargin_3->size[1] != 0)) {
    numalloc = varargin_3->size[0];
    ncolk = varargin_3->size[1];
    for (col = 0; col < ncolk; col++) {
      for (row = 0; row < numalloc; row++) {
        dk = varargin_3->data[row + varargin_3->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cnvardim - 1] = row + 1;
          c->d->data[cnvardim - 1] = dk;
          cnvardim++;
        }
      }

      ccolidx++;
      c->colidx->data[ccolidx] = cnvardim;
    }
  }
}

void g_sparse_horzcat(const emxArray_real_T *varargin_1, int varargin_2_m, int
                      varargin_2_n, const emxArray_real_T *varargin_3, d_sparse *
                      c)
{
  double dk;
  int cidx;
  int cnvardim;
  int col;
  int nrowk;
  int numalloc;
  int row;
  boolean_T allEmpty;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1->size[0];
  isAcceptableEmpty = ((varargin_1->size[0] == 0) || (varargin_1->size[1] == 0));
  foundSize = !isAcceptableEmpty;
  if ((varargin_2_m == 0) || (varargin_2_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  allEmpty = (isAcceptableEmpty && isAcceptableEmpty_tmp);
  if ((!isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_2_m;
  }

  isAcceptableEmpty = ((varargin_3->size[0] == 0) || (varargin_3->size[1] == 0));
  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    c->m = varargin_3->size[0];
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0))) {
    numalloc = 0;
    col = varargin_1->size[0] * varargin_1->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_1->data[cidx] != 0.0) {
        numalloc++;
      }
    }

    cnvardim = varargin_1->size[1];
  }

  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim += varargin_2_n;
  }

  if (allEmpty || ((varargin_3->size[0] != 0) && (varargin_3->size[1] != 0))) {
    nrowk = 0;
    col = varargin_3->size[0] * varargin_3->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_3->data[cidx] != 0.0) {
        nrowk++;
      }
    }

    numalloc += nrowk;
    cnvardim += varargin_3->size[1];
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  col = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, col);
  for (col = 0; col < numalloc; col++) {
    c->d->data[col] = 0.0;
  }

  c->maxnz = numalloc;
  col = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, col);
  c->colidx->data[0] = 1;
  col = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, col);
  for (col = 0; col < numalloc; col++) {
    c->rowidx->data[col] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cidx = 0;
  numalloc = 0;
  if ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0)) {
    nrowk = varargin_1->size[0];
    cnvardim = varargin_1->size[1];
    cidx = 0;
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        dk = varargin_1->data[row + varargin_1->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx] = row + 1;
          c->d->data[cidx] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx + 1;
    }
  }

  if ((varargin_2_m != 0) && (varargin_2_n != 0)) {
    for (col = 0; col < varargin_2_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx + 1;
    }
  }

  if ((varargin_3->size[0] != 0) && (varargin_3->size[1] != 0)) {
    nrowk = varargin_3->size[0];
    cnvardim = varargin_3->size[1];
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        dk = varargin_3->data[row + varargin_3->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx] = row + 1;
          c->d->data[cidx] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx + 1;
    }
  }
}

void h_sparse_horzcat(int varargin_1_m, int varargin_1_n, const emxArray_real_T *
                      varargin_2, int varargin_3_m, int varargin_3_n, int
                      varargin_4_m, int varargin_4_n, d_sparse *c)
{
  double dk;
  int cidx;
  int cnvardim;
  int col;
  int ncolk;
  int numalloc;
  int row;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T c_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1_m;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  foundSize = !isAcceptableEmpty_tmp;
  isAcceptableEmpty = ((varargin_2->size[0] == 0) || (varargin_2->size[1] == 0));
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_2->size[0];
  }

  if ((varargin_3_m == 0) || (varargin_3_n == 0)) {
    b_isAcceptableEmpty_tmp = true;
  } else {
    b_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      c->m = varargin_3_m;
    }
  }

  if ((varargin_4_m == 0) || (varargin_4_n == 0)) {
    c_isAcceptableEmpty_tmp = true;
  } else {
    c_isAcceptableEmpty_tmp = false;
  }

  isAcceptableEmpty = (isAcceptableEmpty_tmp && isAcceptableEmpty &&
                       b_isAcceptableEmpty_tmp && c_isAcceptableEmpty_tmp);
  if ((!c_isAcceptableEmpty_tmp) && (!foundSize)) {
    c->m = varargin_4_m;
  }

  numalloc = 0;
  cnvardim = 0;
  if (isAcceptableEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim = varargin_1_n;
  }

  if (isAcceptableEmpty || ((varargin_2->size[0] != 0) && (varargin_2->size[1]
        != 0))) {
    numalloc = 0;
    cidx = varargin_2->size[0] * varargin_2->size[1];
    for (ncolk = 0; ncolk < cidx; ncolk++) {
      if (varargin_2->data[ncolk] != 0.0) {
        numalloc++;
      }
    }

    cnvardim += varargin_2->size[1];
  }

  if (isAcceptableEmpty || (!b_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_3_n;
  }

  if (isAcceptableEmpty || (!c_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_4_n;
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  cidx = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, cidx);
  for (cidx = 0; cidx < numalloc; cidx++) {
    c->d->data[cidx] = 0.0;
  }

  c->maxnz = numalloc;
  cidx = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, cidx);
  c->colidx->data[0] = 1;
  cidx = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, cidx);
  for (cidx = 0; cidx < numalloc; cidx++) {
    c->rowidx->data[cidx] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cidx = 1;
  numalloc = 0;
  if ((varargin_1_m != 0) && (varargin_1_n != 0)) {
    for (col = 0; col < varargin_1_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = 1;
    }
  }

  if ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0)) {
    cnvardim = varargin_2->size[0];
    ncolk = varargin_2->size[1];
    cidx = 1;
    for (col = 0; col < ncolk; col++) {
      for (row = 0; row < cnvardim; row++) {
        dk = varargin_2->data[row + varargin_2->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_3_m != 0) && (varargin_3_n != 0)) {
    for (col = 0; col < varargin_3_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_4_m != 0) && (varargin_4_n != 0)) {
    for (col = 0; col < varargin_4_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }
}

void i_sparse_horzcat(int varargin_1_m, int varargin_1_n, const emxArray_real_T *
                      varargin_2, int varargin_3_m, int varargin_3_n, const
                      emxArray_real_T *varargin_4, int varargin_5_m, int
                      varargin_5_n, int varargin_6_m, int varargin_6_n, d_sparse
                      *c)
{
  int cidx;
  int cnvardim;
  int col;
  int nrowk;
  int numalloc;
  int row;
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T c_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1_m;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  foundSize = !isAcceptableEmpty_tmp;
  isAcceptableEmpty = ((varargin_2->size[0] == 0) || (varargin_2->size[1] == 0));
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_2->size[0];
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
    c->m = varargin_3_m;
  }

  isAcceptableEmpty = ((varargin_4->size[0] == 0) || (varargin_4->size[1] == 0));
  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_4->size[0];
  }

  if ((varargin_5_m == 0) || (varargin_5_n == 0)) {
    isAcceptableEmpty = true;
  } else {
    isAcceptableEmpty = false;
  }

  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_5_m;
  }

  if ((varargin_6_m == 0) || (varargin_6_n == 0)) {
    c_isAcceptableEmpty_tmp = true;
  } else {
    c_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && c_isAcceptableEmpty_tmp);
  if ((!c_isAcceptableEmpty_tmp) && (!foundSize)) {
    c->m = varargin_6_m;
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim = varargin_1_n;
  }

  if (allEmpty || ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0))) {
    numalloc = 0;
    col = varargin_2->size[0] * varargin_2->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_2->data[cidx] != 0.0) {
        numalloc++;
      }
    }

    cnvardim += varargin_2->size[1];
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_3_n;
  }

  if (allEmpty || ((varargin_4->size[0] != 0) && (varargin_4->size[1] != 0))) {
    nrowk = 0;
    col = varargin_4->size[0] * varargin_4->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_4->data[cidx] != 0.0) {
        nrowk++;
      }
    }

    numalloc += nrowk;
    cnvardim += varargin_4->size[1];
  }

  if (allEmpty || (!isAcceptableEmpty)) {
    cnvardim += varargin_5_n;
  }

  if (allEmpty || (!c_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_6_n;
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  col = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, col);
  for (col = 0; col < numalloc; col++) {
    c->d->data[col] = 0.0;
  }

  c->maxnz = numalloc;
  col = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, col);
  c->colidx->data[0] = 1;
  col = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, col);
  for (col = 0; col < numalloc; col++) {
    c->rowidx->data[col] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cidx = 1;
  numalloc = 0;
  if ((varargin_1_m != 0) && (varargin_1_n != 0)) {
    for (col = 0; col < varargin_1_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = 1;
    }
  }

  if ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0)) {
    nrowk = varargin_2->size[0];
    cnvardim = varargin_2->size[1];
    cidx = 1;
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        if (varargin_2->data[row + varargin_2->size[0] * col] != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = 1.0;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_3_m != 0) && (varargin_3_n != 0)) {
    for (col = 0; col < varargin_3_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_4->size[0] != 0) && (varargin_4->size[1] != 0)) {
    nrowk = varargin_4->size[0];
    cnvardim = varargin_4->size[1];
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        if (varargin_4->data[row + varargin_4->size[0] * col] != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = -1.0;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_5_m != 0) && (varargin_5_n != 0)) {
    for (col = 0; col < varargin_5_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_6_m != 0) && (varargin_6_n != 0)) {
    for (col = 0; col < varargin_6_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }
}

void j_sparse_horzcat(int varargin_1_m, int varargin_1_n, const emxArray_real_T *
                      varargin_2, const emxArray_real_T *varargin_3, int
                      varargin_4_m, int varargin_4_n, int varargin_5_m, int
                      varargin_5_n, d_sparse *c)
{
  double dk;
  int cidx;
  int cnvardim;
  int col;
  int nrowk;
  int numalloc;
  int row;
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1_m;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  foundSize = !isAcceptableEmpty_tmp;
  isAcceptableEmpty = ((varargin_2->size[0] == 0) || (varargin_2->size[1] == 0));
  allEmpty = (isAcceptableEmpty_tmp && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_2->size[0];
  }

  isAcceptableEmpty = ((varargin_3->size[0] == 0) || (varargin_3->size[1] == 0));
  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_3->size[0];
  }

  if ((varargin_4_m == 0) || (varargin_4_n == 0)) {
    isAcceptableEmpty = true;
  } else {
    isAcceptableEmpty = false;
  }

  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_4_m;
  }

  if ((varargin_5_m == 0) || (varargin_5_n == 0)) {
    b_isAcceptableEmpty_tmp = true;
  } else {
    b_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && b_isAcceptableEmpty_tmp);
  if ((!b_isAcceptableEmpty_tmp) && (!foundSize)) {
    c->m = varargin_5_m;
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim = varargin_1_n;
  }

  if (allEmpty || ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0))) {
    numalloc = 0;
    col = varargin_2->size[0] * varargin_2->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_2->data[cidx] != 0.0) {
        numalloc++;
      }
    }

    cnvardim += varargin_2->size[1];
  }

  if (allEmpty || ((varargin_3->size[0] != 0) && (varargin_3->size[1] != 0))) {
    nrowk = 0;
    col = varargin_3->size[0] * varargin_3->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_3->data[cidx] != 0.0) {
        nrowk++;
      }
    }

    numalloc += nrowk;
    cnvardim += varargin_3->size[1];
  }

  if (allEmpty || (!isAcceptableEmpty)) {
    cnvardim += varargin_4_n;
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_5_n;
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  col = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, col);
  for (col = 0; col < numalloc; col++) {
    c->d->data[col] = 0.0;
  }

  c->maxnz = numalloc;
  col = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, col);
  c->colidx->data[0] = 1;
  col = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, col);
  for (col = 0; col < numalloc; col++) {
    c->rowidx->data[col] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cidx = 1;
  numalloc = 0;
  if ((varargin_1_m != 0) && (varargin_1_n != 0)) {
    for (col = 0; col < varargin_1_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = 1;
    }
  }

  if ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0)) {
    nrowk = varargin_2->size[0];
    cnvardim = varargin_2->size[1];
    cidx = 1;
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        dk = varargin_2->data[row + varargin_2->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_3->size[0] != 0) && (varargin_3->size[1] != 0)) {
    nrowk = varargin_3->size[0];
    cnvardim = varargin_3->size[1];
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        dk = varargin_3->data[row + varargin_3->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_4_m != 0) && (varargin_4_n != 0)) {
    for (col = 0; col < varargin_4_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_5_m != 0) && (varargin_5_n != 0)) {
    for (col = 0; col < varargin_5_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }
}

void k_sparse_horzcat(int varargin_1_m, int varargin_1_n, const emxArray_real_T *
                      varargin_2, const emxArray_real_T *varargin_3_d, const
                      emxArray_int32_T *varargin_3_colidx, const
                      emxArray_int32_T *varargin_3_rowidx, int varargin_3_m, int
                      varargin_3_n, int varargin_4_m, int varargin_4_n, int
                      varargin_5_m, int varargin_5_n, d_sparse *c)
{
  double dk;
  int ccolidx;
  int cidx;
  int cnvardim;
  int col;
  int numalloc;
  int row;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T c_isAcceptableEmpty_tmp;
  boolean_T d_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1_m;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  foundSize = !isAcceptableEmpty_tmp;
  isAcceptableEmpty = ((varargin_2->size[0] == 0) || (varargin_2->size[1] == 0));
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_2->size[0];
  }

  if ((varargin_3_m == 0) || (varargin_3_n == 0)) {
    b_isAcceptableEmpty_tmp = true;
  } else {
    b_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      c->m = varargin_3_m;
    }
  }

  if ((varargin_4_m == 0) || (varargin_4_n == 0)) {
    c_isAcceptableEmpty_tmp = true;
  } else {
    c_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      c->m = varargin_4_m;
    }
  }

  if ((varargin_5_m == 0) || (varargin_5_n == 0)) {
    d_isAcceptableEmpty_tmp = true;
  } else {
    d_isAcceptableEmpty_tmp = false;
  }

  isAcceptableEmpty = (isAcceptableEmpty_tmp && isAcceptableEmpty &&
                       b_isAcceptableEmpty_tmp && c_isAcceptableEmpty_tmp &&
                       d_isAcceptableEmpty_tmp);
  if ((!d_isAcceptableEmpty_tmp) && (!foundSize)) {
    c->m = varargin_5_m;
  }

  numalloc = 0;
  cnvardim = 0;
  if (isAcceptableEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim = varargin_1_n;
  }

  if (isAcceptableEmpty || ((varargin_2->size[0] != 0) && (varargin_2->size[1]
        != 0))) {
    numalloc = 0;
    cidx = varargin_2->size[0] * varargin_2->size[1];
    for (ccolidx = 0; ccolidx < cidx; ccolidx++) {
      if (varargin_2->data[ccolidx] != 0.0) {
        numalloc++;
      }
    }

    cnvardim += varargin_2->size[1];
  }

  if (isAcceptableEmpty || (!b_isAcceptableEmpty_tmp)) {
    numalloc = (numalloc + varargin_3_colidx->data[varargin_3_colidx->size[0] -
                1]) - 1;
    cnvardim += varargin_3_n;
  }

  if (isAcceptableEmpty || (!c_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_4_n;
  }

  if (isAcceptableEmpty || (!d_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_5_n;
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  cidx = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, cidx);
  for (cidx = 0; cidx < numalloc; cidx++) {
    c->d->data[cidx] = 0.0;
  }

  c->maxnz = numalloc;
  cidx = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, cidx);
  c->colidx->data[0] = 1;
  cidx = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, cidx);
  for (cidx = 0; cidx < numalloc; cidx++) {
    c->rowidx->data[cidx] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cnvardim = 0;
  ccolidx = 0;
  if ((varargin_1_m != 0) && (varargin_1_n != 0)) {
    for (col = 0; col < varargin_1_n; col++) {
      ccolidx++;
      c->colidx->data[ccolidx] = 1;
    }
  }

  if ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0)) {
    numalloc = varargin_2->size[0];
    cnvardim = varargin_2->size[1];
    cidx = 0;
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < numalloc; row++) {
        dk = varargin_2->data[row + varargin_2->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx] = row + 1;
          c->d->data[cidx] = dk;
          cidx++;
        }
      }

      ccolidx++;
      c->colidx->data[ccolidx] = cidx + 1;
    }

    cnvardim = cidx;
  }

  if ((varargin_3_m != 0) && (varargin_3_n != 0)) {
    cidx = cnvardim - 1;
    numalloc = varargin_3_colidx->data[varargin_3_colidx->size[0] - 1];
    for (row = 0; row <= numalloc - 2; row++) {
      cidx++;
      c->rowidx->data[cidx] = varargin_3_rowidx->data[row];
      c->d->data[cidx] = varargin_3_d->data[row];
    }

    for (col = 0; col < varargin_3_n; col++) {
      ccolidx++;
      c->colidx->data[ccolidx] = varargin_3_colidx->data[col + 1] + cnvardim;
    }

    cnvardim = (cnvardim + varargin_3_colidx->data[varargin_3_colidx->size[0] -
                1]) - 1;
  }

  if ((varargin_4_m != 0) && (varargin_4_n != 0)) {
    for (col = 0; col < varargin_4_n; col++) {
      ccolidx++;
      c->colidx->data[ccolidx] = cnvardim + 1;
    }
  }

  if ((varargin_5_m != 0) && (varargin_5_n != 0)) {
    for (col = 0; col < varargin_5_n; col++) {
      ccolidx++;
      c->colidx->data[ccolidx] = cnvardim + 1;
    }
  }
}

void l_sparse_horzcat(int varargin_1_m, int varargin_1_n, const emxArray_real_T *
                      varargin_2_d, const emxArray_int32_T *varargin_2_colidx,
                      const emxArray_int32_T *varargin_2_rowidx, int
                      varargin_2_m, int varargin_2_n, const emxArray_real_T
                      *varargin_3_d, const emxArray_int32_T *varargin_3_colidx,
                      const emxArray_int32_T *varargin_3_rowidx, int
                      varargin_3_m, int varargin_3_n, int varargin_4_m, int
                      varargin_4_n, d_sparse *c)
{
  int ccolidx;
  int cidx;
  int cnvardim;
  int idx;
  int numalloc;
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T c_isAcceptableEmpty_tmp;
  boolean_T d_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1_m;
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
      c->m = varargin_2_m;
    }
  }

  if ((varargin_3_m == 0) || (varargin_3_n == 0)) {
    c_isAcceptableEmpty_tmp = true;
  } else {
    c_isAcceptableEmpty_tmp = false;
    if (!foundSize) {
      foundSize = true;
      c->m = varargin_3_m;
    }
  }

  if ((varargin_4_m == 0) || (varargin_4_n == 0)) {
    d_isAcceptableEmpty_tmp = true;
  } else {
    d_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (isAcceptableEmpty_tmp && b_isAcceptableEmpty_tmp &&
              c_isAcceptableEmpty_tmp && d_isAcceptableEmpty_tmp);
  if ((!d_isAcceptableEmpty_tmp) && (!foundSize)) {
    c->m = varargin_4_m;
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim = varargin_1_n;
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    numalloc = varargin_2_colidx->data[varargin_2_colidx->size[0] - 1] - 1;
    cnvardim += varargin_2_n;
  }

  if (allEmpty || (!c_isAcceptableEmpty_tmp)) {
    numalloc = (numalloc + varargin_3_colidx->data[varargin_3_colidx->size[0] -
                1]) - 1;
    cnvardim += varargin_3_n;
  }

  if (allEmpty || (!d_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_4_n;
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  cidx = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, cidx);
  for (cidx = 0; cidx < numalloc; cidx++) {
    c->d->data[cidx] = 0.0;
  }

  c->maxnz = numalloc;
  cidx = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, cidx);
  c->colidx->data[0] = 1;
  cidx = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, cidx);
  for (cidx = 0; cidx < numalloc; cidx++) {
    c->rowidx->data[cidx] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cnvardim = 1;
  ccolidx = 0;
  if ((varargin_1_m != 0) && (varargin_1_n != 0)) {
    for (numalloc = 0; numalloc < varargin_1_n; numalloc++) {
      ccolidx++;
      c->colidx->data[ccolidx] = 1;
    }
  }

  if ((varargin_2_m != 0) && (varargin_2_n != 0)) {
    cidx = -1;
    numalloc = varargin_2_colidx->data[varargin_2_colidx->size[0] - 1];
    for (idx = 0; idx <= numalloc - 2; idx++) {
      cidx++;
      c->rowidx->data[cidx] = varargin_2_rowidx->data[idx];
      c->d->data[cidx] = varargin_2_d->data[idx];
    }

    for (numalloc = 0; numalloc < varargin_2_n; numalloc++) {
      ccolidx++;
      c->colidx->data[ccolidx] = varargin_2_colidx->data[numalloc + 1];
    }

    cnvardim = varargin_2_colidx->data[varargin_2_colidx->size[0] - 1];
  }

  if ((varargin_3_m != 0) && (varargin_3_n != 0)) {
    cidx = cnvardim - 2;
    numalloc = varargin_3_colidx->data[varargin_3_colidx->size[0] - 1];
    for (idx = 0; idx <= numalloc - 2; idx++) {
      cidx++;
      c->rowidx->data[cidx] = varargin_3_rowidx->data[idx];
      c->d->data[cidx] = varargin_3_d->data[idx];
    }

    for (numalloc = 0; numalloc < varargin_3_n; numalloc++) {
      ccolidx++;
      c->colidx->data[ccolidx] = (varargin_3_colidx->data[numalloc + 1] +
        cnvardim) - 1;
    }

    cnvardim = (cnvardim + varargin_3_colidx->data[varargin_3_colidx->size[0] -
                1]) - 1;
  }

  if ((varargin_4_m != 0) && (varargin_4_n != 0)) {
    for (numalloc = 0; numalloc < varargin_4_n; numalloc++) {
      ccolidx++;
      c->colidx->data[ccolidx] = cnvardim;
    }
  }
}

void m_sparse_horzcat(const emxArray_real_T *varargin_1, int varargin_2_m, int
                      varargin_2_n, d_sparse *c)
{
  int cidx;
  int cnvardim;
  int col;
  int ncolk;
  int numalloc;
  int row;
  boolean_T allEmpty;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1->size[0];
  isAcceptableEmpty = ((varargin_1->size[0] == 0) || (varargin_1->size[1] == 0));
  if ((varargin_2_m == 0) || (varargin_2_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  allEmpty = (isAcceptableEmpty && isAcceptableEmpty_tmp);
  if ((!isAcceptableEmpty_tmp) && isAcceptableEmpty) {
    c->m = varargin_2_m;
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0))) {
    numalloc = 0;
    ncolk = varargin_1->size[0] * varargin_1->size[1];
    for (cnvardim = 0; cnvardim < ncolk; cnvardim++) {
      if (varargin_1->data[cnvardim] != 0.0) {
        numalloc++;
      }
    }

    cnvardim = varargin_1->size[1];
  }

  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim += varargin_2_n;
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  ncolk = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, ncolk);
  for (ncolk = 0; ncolk < numalloc; ncolk++) {
    c->d->data[ncolk] = 0.0;
  }

  c->maxnz = numalloc;
  ncolk = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, ncolk);
  c->colidx->data[0] = 1;
  ncolk = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, ncolk);
  for (ncolk = 0; ncolk < numalloc; ncolk++) {
    c->rowidx->data[ncolk] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cidx = 0;
  numalloc = 0;
  if ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0)) {
    cnvardim = varargin_1->size[0];
    ncolk = varargin_1->size[1];
    cidx = 0;
    for (col = 0; col < ncolk; col++) {
      for (row = 0; row < cnvardim; row++) {
        if (varargin_1->data[row + varargin_1->size[0] * col] != 0.0) {
          c->rowidx->data[cidx] = row + 1;
          c->d->data[cidx] = 1.0;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx + 1;
    }
  }

  if ((varargin_2_m != 0) && (varargin_2_n != 0)) {
    for (col = 0; col < varargin_2_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx + 1;
    }
  }
}

void n_sparse_horzcat(int varargin_1_m, int varargin_1_n, const emxArray_real_T *
                      varargin_2, const emxArray_real_T *varargin_3, int
                      varargin_4_m, int varargin_4_n, d_sparse *c)
{
  double dk;
  int cidx;
  int cnvardim;
  int col;
  int nrowk;
  int numalloc;
  int row;
  boolean_T allEmpty;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1_m;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  foundSize = !isAcceptableEmpty_tmp;
  isAcceptableEmpty = ((varargin_2->size[0] == 0) || (varargin_2->size[1] == 0));
  allEmpty = (isAcceptableEmpty_tmp && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_2->size[0];
  }

  isAcceptableEmpty = ((varargin_3->size[0] == 0) || (varargin_3->size[1] == 0));
  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_3->size[0];
  }

  if ((varargin_4_m == 0) || (varargin_4_n == 0)) {
    isAcceptableEmpty = true;
  } else {
    isAcceptableEmpty = false;
  }

  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    c->m = varargin_4_m;
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim = varargin_1_n;
  }

  if (allEmpty || ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0))) {
    numalloc = 0;
    col = varargin_2->size[0] * varargin_2->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_2->data[cidx] != 0.0) {
        numalloc++;
      }
    }

    cnvardim += varargin_2->size[1];
  }

  if (allEmpty || ((varargin_3->size[0] != 0) && (varargin_3->size[1] != 0))) {
    nrowk = 0;
    col = varargin_3->size[0] * varargin_3->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_3->data[cidx] != 0.0) {
        nrowk++;
      }
    }

    numalloc += nrowk;
    cnvardim += varargin_3->size[1];
  }

  if (allEmpty || (!isAcceptableEmpty)) {
    cnvardim += varargin_4_n;
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  col = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, col);
  for (col = 0; col < numalloc; col++) {
    c->d->data[col] = 0.0;
  }

  c->maxnz = numalloc;
  col = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, col);
  c->colidx->data[0] = 1;
  col = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, col);
  for (col = 0; col < numalloc; col++) {
    c->rowidx->data[col] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cidx = 1;
  numalloc = 0;
  if ((varargin_1_m != 0) && (varargin_1_n != 0)) {
    for (col = 0; col < varargin_1_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = 1;
    }
  }

  if ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0)) {
    nrowk = varargin_2->size[0];
    cnvardim = varargin_2->size[1];
    cidx = 1;
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        dk = varargin_2->data[row + varargin_2->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_3->size[0] != 0) && (varargin_3->size[1] != 0)) {
    nrowk = varargin_3->size[0];
    cnvardim = varargin_3->size[1];
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        dk = varargin_3->data[row + varargin_3->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_4_m != 0) && (varargin_4_n != 0)) {
    for (col = 0; col < varargin_4_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }
}

void o_sparse_horzcat(int K,int varargin_1_m, int varargin_1_n, const emxArray_real_T *
                      varargin_2, int varargin_3_m, int varargin_3_n, d_sparse
                      *c)
{
  double dk;
  int cidx;
  int cnvardim;
  int col;
  int numalloc;
  int row;
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1_m;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  foundSize = !isAcceptableEmpty_tmp;
  allEmpty = (varargin_2->size[1] == 0);
  if ((!allEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = K;
  }

  if ((varargin_3_m == 0) || (varargin_3_n == 0)) {
    b_isAcceptableEmpty_tmp = true;
  } else {
    b_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (isAcceptableEmpty_tmp && allEmpty && b_isAcceptableEmpty_tmp);
  if ((!b_isAcceptableEmpty_tmp) && (!foundSize)) {
    c->m = varargin_3_m;
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim = varargin_1_n;
  }

  if (allEmpty || (varargin_2->size[1] != 0)) {
    numalloc = 0;
    col = K * varargin_2->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_2->data[cidx] != 0.0) {
        numalloc++;
      }
    }

    cnvardim += varargin_2->size[1];
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_3_n;
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  col = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, col);
  for (col = 0; col < numalloc; col++) {
    c->d->data[col] = 0.0;
  }

  c->maxnz = numalloc;
  col = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, col);
  c->colidx->data[0] = 1;
  col = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, col);
  for (col = 0; col < numalloc; col++) {
    c->rowidx->data[col] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cidx = 0;
  numalloc = 0;
  if ((varargin_1_m != 0) && (varargin_1_n != 0)) {
    for (col = 0; col < varargin_1_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = 1;
    }
  }

  if (varargin_2->size[1] != 0) {
    cnvardim = varargin_2->size[1];
    cidx = 0;
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < K; row++) {
        dk = varargin_2->data[row + K * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx] = row + 1;
          c->d->data[cidx] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx + 1;
    }
  }

  if ((varargin_3_m != 0) && (varargin_3_n != 0)) {
    for (col = 0; col < varargin_3_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx + 1;
    }
  }
}

void p_sparse_horzcat(int K,int varargin_1_m, int varargin_1_n, const emxArray_real_T *
                      varargin_2_d, const emxArray_int32_T *varargin_2_colidx,
                      const emxArray_int32_T *varargin_2_rowidx, int
                      varargin_2_m, int varargin_2_n, const emxArray_real_T
                      *varargin_3, int varargin_4_m, int varargin_4_n, d_sparse *
                      c)
{
  double dk;
  int cnvardim;
  int k;
  int numalloc;
  int nzCount;
  int row;
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T c_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1_m;
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
      c->m = varargin_2_m;
    }
  }

  allEmpty = (varargin_3->size[1] == 0);
  if ((!allEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = K;
  }

  if ((varargin_4_m == 0) || (varargin_4_n == 0)) {
    c_isAcceptableEmpty_tmp = true;
  } else {
    c_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (isAcceptableEmpty_tmp && b_isAcceptableEmpty_tmp && allEmpty &&
              c_isAcceptableEmpty_tmp);
  if ((!c_isAcceptableEmpty_tmp) && (!foundSize)) {
    c->m = varargin_4_m;
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim = varargin_1_n;
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    numalloc = varargin_2_colidx->data[varargin_2_colidx->size[0] - 1] - 1;
    cnvardim += varargin_2_n;
  }

  if (allEmpty || (varargin_3->size[1] != 0)) {
    nzCount = 0;
    row = K * varargin_3->size[1];
    for (k = 0; k < row; k++) {
      if (varargin_3->data[k] != 0.0) {
        nzCount++;
      }
    }

    numalloc += nzCount;
    cnvardim += varargin_3->size[1];
  }

  if (allEmpty || (!c_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_4_n;
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  row = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, row);
  for (row = 0; row < numalloc; row++) {
    c->d->data[row] = 0.0;
  }

  c->maxnz = numalloc;
  row = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, row);
  c->colidx->data[0] = 1;
  row = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, row);
  for (row = 0; row < numalloc; row++) {
    c->rowidx->data[row] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  nzCount = 1;
  cnvardim = 0;
  if ((varargin_1_m != 0) && (varargin_1_n != 0)) {
    for (k = 0; k < varargin_1_n; k++) {
      cnvardim++;
      c->colidx->data[cnvardim] = 1;
    }
  }

  if ((varargin_2_m != 0) && (varargin_2_n != 0)) {
    nzCount = -1;
    numalloc = varargin_2_colidx->data[varargin_2_colidx->size[0] - 1];
    for (k = 0; k <= numalloc - 2; k++) {
      nzCount++;
      c->rowidx->data[nzCount] = varargin_2_rowidx->data[k];
      c->d->data[nzCount] = varargin_2_d->data[k];
    }

    for (k = 0; k < varargin_2_n; k++) {
      cnvardim++;
      c->colidx->data[cnvardim] = varargin_2_colidx->data[k + 1];
    }

    nzCount = varargin_2_colidx->data[varargin_2_colidx->size[0] - 1];
  }

  if (varargin_3->size[1] != 0) {
    numalloc = varargin_3->size[1];
    for (k = 0; k < numalloc; k++) {
      for (row = 0; row < K; row++) {
        dk = varargin_3->data[row + K * k];
        if (dk != 0.0) {
          c->rowidx->data[nzCount - 1] = row + 1;
          c->d->data[nzCount - 1] = dk;
          nzCount++;
        }
      }

      cnvardim++;
      c->colidx->data[cnvardim] = nzCount;
    }
  }

  if ((varargin_4_m != 0) && (varargin_4_n != 0)) {
    for (k = 0; k < varargin_4_n; k++) {
      cnvardim++;
      c->colidx->data[cnvardim] = nzCount;
    }
  }
}

void q_sparse_horzcat(int varargin_1_m, int varargin_1_n, const emxArray_real_T *
                      varargin_2, int varargin_3_m, int varargin_3_n, const
                      emxArray_real_T *varargin_4, d_sparse *c)
{
  double dk;
  int cidx;
  int cnvardim;
  int col;
  int nrowk;
  int numalloc;
  int row;
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1_m;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  foundSize = !isAcceptableEmpty_tmp;
  isAcceptableEmpty = ((varargin_2->size[0] == 0) || (varargin_2->size[1] == 0));
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_2->size[0];
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
    c->m = varargin_3_m;
  }

  isAcceptableEmpty = ((varargin_4->size[0] == 0) || (varargin_4->size[1] == 0));
  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    c->m = varargin_4->size[0];
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim = varargin_1_n;
  }

  if (allEmpty || ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0))) {
    numalloc = 0;
    col = varargin_2->size[0] * varargin_2->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_2->data[cidx] != 0.0) {
        numalloc++;
      }
    }

    cnvardim += varargin_2->size[1];
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_3_n;
  }

  if (allEmpty || ((varargin_4->size[0] != 0) && (varargin_4->size[1] != 0))) {
    nrowk = 0;
    col = varargin_4->size[0] * varargin_4->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_4->data[cidx] != 0.0) {
        nrowk++;
      }
    }

    numalloc += nrowk;
    cnvardim += varargin_4->size[1];
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  col = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, col);
  for (col = 0; col < numalloc; col++) {
    c->d->data[col] = 0.0;
  }

  c->maxnz = numalloc;
  col = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, col);
  c->colidx->data[0] = 1;
  col = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, col);
  for (col = 0; col < numalloc; col++) {
    c->rowidx->data[col] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cidx = 0;
  numalloc = 0;
  if ((varargin_1_m != 0) && (varargin_1_n != 0)) {
    for (col = 0; col < varargin_1_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = 1;
    }
  }

  if ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0)) {
    nrowk = varargin_2->size[0];
    cnvardim = varargin_2->size[1];
    cidx = 0;
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        dk = varargin_2->data[row + varargin_2->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx] = row + 1;
          c->d->data[cidx] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx + 1;
    }
  }

  if ((varargin_3_m != 0) && (varargin_3_n != 0)) {
    for (col = 0; col < varargin_3_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx + 1;
    }
  }

  if ((varargin_4->size[0] != 0) && (varargin_4->size[1] != 0)) {
    nrowk = varargin_4->size[0];
    cnvardim = varargin_4->size[1];
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        dk = varargin_4->data[row + varargin_4->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx] = row + 1;
          c->d->data[cidx] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx + 1;
    }
  }
}

void r_sparse_horzcat(const emxArray_real_T *varargin_1_d, const
                      emxArray_int32_T *varargin_1_colidx, const
                      emxArray_int32_T *varargin_1_rowidx, int varargin_1_m, int
                      varargin_1_n, int varargin_2_m, int varargin_2_n, d_sparse
                      *c)
{
  int ccolidx;
  int cnvardim;
  int idx;
  int numalloc;
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1_m;
  if ((varargin_1_m == 0) || (varargin_1_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  if ((varargin_2_m == 0) || (varargin_2_n == 0)) {
    b_isAcceptableEmpty_tmp = true;
  } else {
    b_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (isAcceptableEmpty_tmp && b_isAcceptableEmpty_tmp);
  if ((!b_isAcceptableEmpty_tmp) && isAcceptableEmpty_tmp) {
    c->m = varargin_2_m;
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    numalloc = varargin_1_colidx->data[varargin_1_colidx->size[0] - 1] - 1;
    cnvardim = varargin_1_n;
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_2_n;
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  ccolidx = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, ccolidx);
  for (ccolidx = 0; ccolidx < numalloc; ccolidx++) {
    c->d->data[ccolidx] = 0.0;
  }

  c->maxnz = numalloc;
  ccolidx = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, ccolidx);
  c->colidx->data[0] = 1;
  ccolidx = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, ccolidx);
  for (ccolidx = 0; ccolidx < numalloc; ccolidx++) {
    c->rowidx->data[ccolidx] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  numalloc = 1;
  ccolidx = 0;
  if ((varargin_1_m != 0) && (varargin_1_n != 0)) {
    numalloc = -1;
    cnvardim = varargin_1_colidx->data[varargin_1_colidx->size[0] - 1];
    for (idx = 0; idx <= cnvardim - 2; idx++) {
      numalloc++;
      c->rowidx->data[numalloc] = varargin_1_rowidx->data[idx];
      c->d->data[numalloc] = varargin_1_d->data[idx];
    }

    for (cnvardim = 0; cnvardim < varargin_1_n; cnvardim++) {
      ccolidx++;
      c->colidx->data[ccolidx] = varargin_1_colidx->data[cnvardim + 1];
    }

    numalloc = varargin_1_colidx->data[varargin_1_colidx->size[0] - 1];
  }

  if ((varargin_2_m != 0) && (varargin_2_n != 0)) {
    for (cnvardim = 0; cnvardim < varargin_2_n; cnvardim++) {
      ccolidx++;
      c->colidx->data[ccolidx] = numalloc;
    }
  }
}

void sparse_horzcat(const emxArray_real_T *varargin_1, int varargin_2_m, int
                    varargin_2_n, const emxArray_real_T *varargin_3, int
                    varargin_4_m, int varargin_4_n, int varargin_5_m, int
                    varargin_5_n, int varargin_6_m, int varargin_6_n, const
                    emxArray_real_T *varargin_7, int varargin_8_m, int
                    varargin_8_n, int varargin_9_m, int varargin_9_n, d_sparse
                    *c)
{
  double dk;
  int cidx;
  int cnvardim;
  int col;
  int nrowk;
  int numalloc;
  int row;
  boolean_T allEmpty;
  boolean_T b_isAcceptableEmpty_tmp;
  boolean_T c_isAcceptableEmpty_tmp;
  boolean_T d_isAcceptableEmpty_tmp;
  boolean_T e_isAcceptableEmpty_tmp;
  boolean_T foundSize;
  boolean_T isAcceptableEmpty;
  boolean_T isAcceptableEmpty_tmp;
  c->m = varargin_1->size[0];
  isAcceptableEmpty = ((varargin_1->size[0] == 0) || (varargin_1->size[1] == 0));
  foundSize = !isAcceptableEmpty;
  if ((varargin_2_m == 0) || (varargin_2_n == 0)) {
    isAcceptableEmpty_tmp = true;
  } else {
    isAcceptableEmpty_tmp = false;
  }

  allEmpty = (isAcceptableEmpty && isAcceptableEmpty_tmp);
  if ((!isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_2_m;
  }

  isAcceptableEmpty = ((varargin_3->size[0] == 0) || (varargin_3->size[1] == 0));
  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_3->size[0];
  }

  if ((varargin_4_m == 0) || (varargin_4_n == 0)) {
    b_isAcceptableEmpty_tmp = true;
  } else {
    b_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && b_isAcceptableEmpty_tmp);
  if ((!b_isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_4_m;
  }

  if ((varargin_5_m == 0) || (varargin_5_n == 0)) {
    c_isAcceptableEmpty_tmp = true;
  } else {
    c_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && c_isAcceptableEmpty_tmp);
  if ((!c_isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_5_m;
  }

  if ((varargin_6_m == 0) || (varargin_6_n == 0)) {
    d_isAcceptableEmpty_tmp = true;
  } else {
    d_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && d_isAcceptableEmpty_tmp);
  if ((!d_isAcceptableEmpty_tmp) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_6_m;
  }

  isAcceptableEmpty = ((varargin_7->size[0] == 0) || (varargin_7->size[1] == 0));
  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_7->size[0];
  }

  if ((varargin_8_m == 0) || (varargin_8_n == 0)) {
    isAcceptableEmpty = true;
  } else {
    isAcceptableEmpty = false;
  }

  allEmpty = (allEmpty && isAcceptableEmpty);
  if ((!isAcceptableEmpty) && (!foundSize)) {
    foundSize = true;
    c->m = varargin_8_m;
  }

  if ((varargin_9_m == 0) || (varargin_9_n == 0)) {
    e_isAcceptableEmpty_tmp = true;
  } else {
    e_isAcceptableEmpty_tmp = false;
  }

  allEmpty = (allEmpty && e_isAcceptableEmpty_tmp);
  if ((!e_isAcceptableEmpty_tmp) && (!foundSize)) {
    c->m = varargin_9_m;
  }

  numalloc = 0;
  cnvardim = 0;
  if (allEmpty || ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0))) {
    numalloc = 0;
    col = varargin_1->size[0] * varargin_1->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_1->data[cidx] != 0.0) {
        numalloc++;
      }
    }

    cnvardim = varargin_1->size[1];
  }

  if (allEmpty || (!isAcceptableEmpty_tmp)) {
    cnvardim += varargin_2_n;
  }

  if (allEmpty || ((varargin_3->size[0] != 0) && (varargin_3->size[1] != 0))) {
    nrowk = 0;
    col = varargin_3->size[0] * varargin_3->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_3->data[cidx] != 0.0) {
        nrowk++;
      }
    }

    numalloc += nrowk;
    cnvardim += varargin_3->size[1];
  }

  if (allEmpty || (!b_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_4_n;
  }

  if (allEmpty || (!c_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_5_n;
  }

  if (allEmpty || (!d_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_6_n;
  }

  if (allEmpty || ((varargin_7->size[0] != 0) && (varargin_7->size[1] != 0))) {
    nrowk = 0;
    col = varargin_7->size[0] * varargin_7->size[1];
    for (cidx = 0; cidx < col; cidx++) {
      if (varargin_7->data[cidx] != 0.0) {
        nrowk++;
      }
    }

    numalloc += nrowk;
    cnvardim += varargin_7->size[1];
  }

  if (allEmpty || (!isAcceptableEmpty)) {
    cnvardim += varargin_8_n;
  }

  if (allEmpty || (!e_isAcceptableEmpty_tmp)) {
    cnvardim += varargin_9_n;
  }

  c->n = cnvardim;
  if (numalloc < 1) {
    numalloc = 1;
  }

  col = c->d->size[0];
  c->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(c->d, col);
  for (col = 0; col < numalloc; col++) {
    c->d->data[col] = 0.0;
  }

  c->maxnz = numalloc;
  col = c->colidx->size[0];
  c->colidx->size[0] = cnvardim + 1;
  emxEnsureCapacity_int32_T(c->colidx, col);
  c->colidx->data[0] = 1;
  col = c->rowidx->size[0];
  c->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(c->rowidx, col);
  for (col = 0; col < numalloc; col++) {
    c->rowidx->data[col] = 0;
  }

  for (numalloc = 0; numalloc < cnvardim; numalloc++) {
    c->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(c);
  cidx = 1;
  numalloc = 0;
  if ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0)) {
    nrowk = varargin_1->size[0];
    cnvardim = varargin_1->size[1];
    cidx = 1;
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        if (varargin_1->data[row + varargin_1->size[0] * col] != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = -1.0;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_2_m != 0) && (varargin_2_n != 0)) {
    for (col = 0; col < varargin_2_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_3->size[0] != 0) && (varargin_3->size[1] != 0)) {
    nrowk = varargin_3->size[0];
    cnvardim = varargin_3->size[1];
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        if (varargin_3->data[row + varargin_3->size[0] * col] != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = -1.0;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_4_m != 0) && (varargin_4_n != 0)) {
    for (col = 0; col < varargin_4_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_5_m != 0) && (varargin_5_n != 0)) {
    for (col = 0; col < varargin_5_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_6_m != 0) && (varargin_6_n != 0)) {
    for (col = 0; col < varargin_6_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_7->size[0] != 0) && (varargin_7->size[1] != 0)) {
    nrowk = varargin_7->size[0];
    cnvardim = varargin_7->size[1];
    for (col = 0; col < cnvardim; col++) {
      for (row = 0; row < nrowk; row++) {
        dk = varargin_7->data[row + varargin_7->size[0] * col];
        if (dk != 0.0) {
          c->rowidx->data[cidx - 1] = row + 1;
          c->d->data[cidx - 1] = dk;
          cidx++;
        }
      }

      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_8_m != 0) && (varargin_8_n != 0)) {
    for (col = 0; col < varargin_8_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }

  if ((varargin_9_m != 0) && (varargin_9_n != 0)) {
    for (col = 0; col < varargin_9_n; col++) {
      numalloc++;
      c->colidx->data[numalloc] = cidx;
    }
  }
}

/* End of code generation (horzcat.c) */
