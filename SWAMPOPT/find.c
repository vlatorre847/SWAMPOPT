/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * find.c
 *
 * Code generation for function 'find'
 *
 */

/* Include files */
#include "find.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_eml_find(const emxArray_int32_T *x_colidx, const emxArray_int32_T
                *x_rowidx, emxArray_int32_T *i, emxArray_int32_T *j)
{
  int col;
  int idx;
  int nx;
  nx = x_colidx->data[x_colidx->size[0] - 1] - 1;
  if (x_colidx->data[x_colidx->size[0] - 1] - 1 == 0) {
    i->size[0] = 0;
    j->size[0] = 0;
  } else {
    col = i->size[0];
    i->size[0] = x_colidx->data[x_colidx->size[0] - 1] - 1;
    emxEnsureCapacity_int32_T(i, col);
    col = j->size[0];
    j->size[0] = x_colidx->data[x_colidx->size[0] - 1] - 1;
    emxEnsureCapacity_int32_T(j, col);
    for (idx = 0; idx < nx; idx++) {
      i->data[idx] = x_rowidx->data[idx];
    }

    idx = 0;
    col = 1;
    while (idx < nx) {
      if (idx == x_colidx->data[col] - 1) {
        col++;
      } else {
        idx++;
        j->data[idx - 1] = col;
      }
    }

    if (x_colidx->data[x_colidx->size[0] - 1] - 1 == 1) {
      if (idx == 0) {
        i->size[0] = 0;
        j->size[0] = 0;
      }
    } else {
      col = i->size[0];
      if (1 > idx) {
        i->size[0] = 0;
      } else {
        i->size[0] = idx;
      }

      emxEnsureCapacity_int32_T(i, col);
      col = j->size[0];
      if (1 > idx) {
        j->size[0] = 0;
      } else {
        j->size[0] = idx;
      }

      emxEnsureCapacity_int32_T(j, col);
    }
  }
}

void eml_find(const emxArray_boolean_T *x, emxArray_int32_T *i)
{
  int idx;
  int ii;
  int nx;
  boolean_T exitg1;
  nx = x->size[0] * x->size[1];
  idx = 0;
  ii = i->size[0];
  i->size[0] = nx;
  emxEnsureCapacity_int32_T(i, ii);
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii <= nx - 1)) {
    if (x->data[ii]) {
      idx++;
      i->data[idx - 1] = ii + 1;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        ii++;
      }
    } else {
      ii++;
    }
  }

  if (nx == 1) {
    if (idx == 0) {
      i->size[0] = 0;
    }
  } else {
    ii = i->size[0];
    if (1 > idx) {
      i->size[0] = 0;
    } else {
      i->size[0] = idx;
    }

    emxEnsureCapacity_int32_T(i, ii);
  }
}

/* End of code generation (find.c) */
