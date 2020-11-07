/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum.c
 *
 * Code generation for function 'sum'
 *
 */

/* Include files */
#include "sum.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_sum(const emxArray_boolean_T *x_d, const emxArray_int32_T *x_colidx, int
           x_m, int x_n, emxArray_real_T *y_d, emxArray_int32_T *y_colidx,
           emxArray_int32_T *y_rowidx, int *y_n)
{
  double r;
  int col;
  int numalloc;
  int sn;
  int xend;
  int xp;
  int xstart;
  if ((x_m == 0) || (x_n == 0) || (x_m == 0)) {
    if (x_n < 0) {
      sn = 0;
    } else {
      sn = x_n;
    }

    xstart = y_colidx->size[0];
    y_colidx->size[0] = sn + 1;
    emxEnsureCapacity_int32_T(y_colidx, xstart);
    for (xstart = 0; xstart <= sn; xstart++) {
      y_colidx->data[xstart] = 1;
    }

    xstart = y_d->size[0];
    y_d->size[0] = 1;
    emxEnsureCapacity_real_T(y_d, xstart);
    y_d->data[0] = 0.0;
    xstart = y_rowidx->size[0];
    y_rowidx->size[0] = 1;
    emxEnsureCapacity_int32_T(y_rowidx, xstart);
    y_rowidx->data[0] = 1;
  } else {
    if (x_n < x_colidx->data[x_colidx->size[0] - 1] - 1) {
      numalloc = x_n;
    } else {
      numalloc = x_colidx->data[x_colidx->size[0] - 1] - 1;
    }

    sn = x_n;
    if (numalloc < 1) {
      numalloc = 1;
    }

    xstart = y_d->size[0];
    y_d->size[0] = numalloc;
    emxEnsureCapacity_real_T(y_d, xstart);
    for (xstart = 0; xstart < numalloc; xstart++) {
      y_d->data[xstart] = 0.0;
    }

    xstart = y_colidx->size[0];
    y_colidx->size[0] = x_n + 1;
    emxEnsureCapacity_int32_T(y_colidx, xstart);
    xstart = y_rowidx->size[0];
    y_rowidx->size[0] = numalloc;
    emxEnsureCapacity_int32_T(y_rowidx, xstart);
    for (xstart = 0; xstart < numalloc; xstart++) {
      y_rowidx->data[xstart] = 0;
    }

    for (numalloc = 0; numalloc < x_n; numalloc++) {
      y_colidx->data[numalloc + 1] = 1;
    }

    xstart = y_colidx->size[0];
    for (numalloc = 0; numalloc <= xstart - 2; numalloc++) {
      y_colidx->data[numalloc] = 1;
    }

    y_colidx->data[y_colidx->size[0] - 1] = 1;
    y_colidx->data[0] = 1;
    numalloc = 1;
    for (col = 0; col < x_n; col++) {
      xstart = x_colidx->data[col];
      xend = x_colidx->data[col + 1] - 1;
      r = 0.0;
      for (xp = xstart; xp <= xend; xp++) {
        r += (double)x_d->data[xp - 1];
      }

      if (r != 0.0) {
        y_d->data[numalloc - 1] = r;
        numalloc++;
      }

      y_colidx->data[col + 1] = numalloc;
    }

    xstart = y_colidx->data[y_colidx->size[0] - 1];
    for (numalloc = 0; numalloc <= xstart - 2; numalloc++) {
      y_rowidx->data[numalloc] = 1;
    }
  }

  *y_n = sn;
}

double sum(const double x_data[], const int x_size[1])
{
  double y;
  int k;
  int vlen;
  vlen = x_size[0];
  if (x_size[0] == 0) {
    y = 0.0;
  } else {
    y = x_data[0];
    for (k = 2; k <= vlen; k++) {
      y += x_data[k - 1];
    }
  }

  return y;
}

/* End of code generation (sum.c) */
