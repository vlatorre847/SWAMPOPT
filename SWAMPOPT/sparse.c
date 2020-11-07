/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sparse.c
 *
 * Code generation for function 'sparse'
 *
 */

/* Include files */
#include "sparse.h"
#include "fillIn.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "introsort.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_sparse(const emxArray_int32_T *varargin_1, const emxArray_int32_T
              *varargin_2, const emxArray_real_T *varargin_3, double varargin_4,
              double varargin_5, d_sparse *y)
{
  cell_wrap_6 this_tunableEnvironment[2];
  emxArray_int32_T *sortedIndices;
  emxArray_int32_T *t;
  int k;
  int nc;
  int ns;
  int ny;
  emxInitMatrix_cell_wrap_6(this_tunableEnvironment);
  nc = varargin_2->size[0];
  ns = varargin_1->size[0];
  k = this_tunableEnvironment[1].f1->size[0];
  this_tunableEnvironment[1].f1->size[0] = varargin_1->size[0];
  emxEnsureCapacity_int32_T(this_tunableEnvironment[1].f1, k);
  for (k = 0; k < ns; k++) {
    this_tunableEnvironment[1].f1->data[k] = varargin_1->data[k];
  }

  ns = varargin_2->size[0];
  k = this_tunableEnvironment[0].f1->size[0];
  this_tunableEnvironment[0].f1->size[0] = varargin_2->size[0];
  emxEnsureCapacity_int32_T(this_tunableEnvironment[0].f1, k);
  for (k = 0; k < ns; k++) {
    this_tunableEnvironment[0].f1->data[k] = varargin_2->data[k];
  }

  emxInit_int32_T(&sortedIndices, 1);
  k = sortedIndices->size[0];
  sortedIndices->size[0] = varargin_2->size[0];
  emxEnsureCapacity_int32_T(sortedIndices, k);
  for (k = 0; k < nc; k++) {
    sortedIndices->data[k] = k + 1;
  }

  emxInit_int32_T(&t, 1);
  introsort(sortedIndices, this_tunableEnvironment[0].f1->size[0],
            this_tunableEnvironment);
  ny = this_tunableEnvironment[0].f1->size[0];
  k = t->size[0];
  t->size[0] = this_tunableEnvironment[0].f1->size[0];
  emxEnsureCapacity_int32_T(t, k);
  ns = this_tunableEnvironment[0].f1->size[0];
  for (k = 0; k < ns; k++) {
    t->data[k] = this_tunableEnvironment[0].f1->data[k];
  }

  for (k = 0; k < ny; k++) {
    this_tunableEnvironment[0].f1->data[k] = t->data[sortedIndices->data[k] - 1];
  }

  ny = this_tunableEnvironment[1].f1->size[0];
  k = t->size[0];
  t->size[0] = this_tunableEnvironment[1].f1->size[0];
  emxEnsureCapacity_int32_T(t, k);
  ns = this_tunableEnvironment[1].f1->size[0];
  for (k = 0; k < ns; k++) {
    t->data[k] = this_tunableEnvironment[1].f1->data[k];
  }

  for (k = 0; k < ny; k++) {
    this_tunableEnvironment[1].f1->data[k] = t->data[sortedIndices->data[k] - 1];
  }

  emxFree_int32_T(&t);
  y->m = (int)varargin_4;
  y->n = (int)varargin_5;
  if (varargin_2->size[0] >= 1) {
    ns = varargin_2->size[0];
  } else {
    ns = 1;
  }

  k = y->d->size[0];
  y->d->size[0] = ns;
  emxEnsureCapacity_real_T(y->d, k);
  for (k = 0; k < ns; k++) {
    y->d->data[k] = 0.0;
  }

  y->maxnz = ns;
  k = y->colidx->size[0];
  y->colidx->size[0] = (int)varargin_5 + 1;
  emxEnsureCapacity_int32_T(y->colidx, k);
  y->colidx->data[0] = 1;
  k = y->rowidx->size[0];
  y->rowidx->size[0] = ns;
  emxEnsureCapacity_int32_T(y->rowidx, k);
  for (k = 0; k < ns; k++) {
    y->rowidx->data[k] = 0;
  }

  ns = 0;
  k = (int)varargin_5;
  for (ny = 0; ny < k; ny++) {
    while ((ns + 1 <= nc) && (this_tunableEnvironment[0].f1->data[ns] == ny + 1))
    {
      y->rowidx->data[ns] = this_tunableEnvironment[1].f1->data[ns];
      ns++;
    }

    y->colidx->data[ny + 1] = ns + 1;
  }

  emxFreeMatrix_cell_wrap_6(this_tunableEnvironment);
  for (k = 0; k < nc; k++) {
    y->d->data[k] = varargin_3->data[sortedIndices->data[k] - 1];
  }

  emxFree_int32_T(&sortedIndices);
  sparse_fillIn(y);
}

void c_sparse(double varargin_2, emxArray_real_T *y_d, emxArray_int32_T
              *y_colidx, emxArray_int32_T *y_rowidx, int *y_n)
{
  int i;
  int loop_ub;
  i = y_d->size[0];
  y_d->size[0] = 1;
  emxEnsureCapacity_real_T(y_d, i);
  y_d->data[0] = 0.0;
  i = y_colidx->size[0];
  y_colidx->size[0] = (int)varargin_2 + 1;
  emxEnsureCapacity_int32_T(y_colidx, i);
  loop_ub = (int)varargin_2;
  for (i = 0; i <= loop_ub; i++) {
    y_colidx->data[i] = 1;
  }

  i = y_rowidx->size[0];
  y_rowidx->size[0] = 1;
  emxEnsureCapacity_int32_T(y_rowidx, i);
  y_rowidx->data[0] = 1;
  *y_n = (int)varargin_2;
}

void sparse(double varargin_1, double varargin_2, emxArray_real_T *y_d,
            emxArray_int32_T *y_colidx, emxArray_int32_T *y_rowidx, int *y_m,
            int *y_n, int *y_maxnz)
{
  int i;
  int loop_ub;
  i = y_d->size[0];
  y_d->size[0] = 1;
  emxEnsureCapacity_real_T(y_d, i);
  y_d->data[0] = 0.0;
  i = y_colidx->size[0];
  y_colidx->size[0] = (int)varargin_2 + 1;
  emxEnsureCapacity_int32_T(y_colidx, i);
  loop_ub = (int)varargin_2;
  for (i = 0; i <= loop_ub; i++) {
    y_colidx->data[i] = 1;
  }

  i = y_rowidx->size[0];
  y_rowidx->size[0] = 1;
  emxEnsureCapacity_int32_T(y_rowidx, i);
  y_rowidx->data[0] = 1;
  *y_m = (int)varargin_1;
  *y_n = (int)varargin_2;
  *y_maxnz = 1;
}

/* End of code generation (sparse.c) */
