/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * diag.c
 *
 * Code generation for function 'diag'
 *
 */

/* Include files */
#include "diag.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void b_diag(int l, const double v[l], double d[l*l])
{
  int j;
  memset(&d[0], 0, l*l * sizeof(double));
  for (j = 0; j < l; j++) {
    d[j + l * j] = v[j];
  }
}

void c_diag(const emxArray_real_T *v, emxArray_real_T *d)
{
  int i;
  int loop_ub;
  int nv;
  nv = v->size[0];
  i = d->size[0] * d->size[1];
  d->size[0] = v->size[0];
  d->size[1] = v->size[0];
  emxEnsureCapacity_real_T(d, i);
  loop_ub = v->size[0] * v->size[0];
  for (i = 0; i < loop_ub; i++) {
    d->data[i] = 0.0;
  }

  for (loop_ub = 0; loop_ub < nv; loop_ub++) {
    d->data[loop_ub + d->size[0] * loop_ub] = v->data[loop_ub];
  }
}

void d_diag(const emxArray_real_T *v, emxArray_real_T *d)
{
  int i;
  int loop_ub;
  int nv;
  nv = v->size[0];
  i = d->size[0] * d->size[1];
  d->size[0] = v->size[0] + 1;
  d->size[1] = v->size[0] + 1;
  emxEnsureCapacity_real_T(d, i);
  loop_ub = (v->size[0] + 1) * (v->size[0] + 1);
  for (i = 0; i < loop_ub; i++) {
    d->data[i] = 0.0;
  }

  for (loop_ub = 0; loop_ub < nv; loop_ub++) {
    d->data[loop_ub + d->size[0] * (loop_ub + 1)] = 1.0;
  }
}

void diag(const emxArray_real_T *v, emxArray_real_T *d)
{
  int i;
  int loop_ub;
  int nv;
  nv = v->size[1];
  i = d->size[0] * d->size[1];
  d->size[0] = v->size[1];
  d->size[1] = v->size[1];
  emxEnsureCapacity_real_T(d, i);
  loop_ub = v->size[1] * v->size[1];
  for (i = 0; i < loop_ub; i++) {
    d->data[i] = 0.0;
  }

  for (loop_ub = 0; loop_ub < nv; loop_ub++) {
    d->data[loop_ub + d->size[0] * loop_ub] = v->data[loop_ub];
  }
}

void e_diag(const emxArray_real_T *v, double K, emxArray_real_T *d)
{
  int absk;
  int i;
  int m;
  int nv;
  boolean_T negk;
  if (K < 0.0) {
    if (-K < 2.147483647E+9) {
      absk = (int)-K;
    } else {
      absk = MAX_int32_T;
    }

    negk = true;
  } else {
    if (K < 2.147483647E+9) {
      absk = (int)K;
    } else {
      absk = MAX_int32_T;
    }

    negk = false;
  }

  nv = v->size[0] - 1;
  m = v->size[0] + absk;
  i = d->size[0] * d->size[1];
  d->size[0] = m;
  d->size[1] = m;
  emxEnsureCapacity_real_T(d, i);
  m *= m;
  for (i = 0; i < m; i++) {
    d->data[i] = 0.0;
  }

  if (negk) {
    for (m = 0; m <= nv; m++) {
      d->data[(m + absk) + d->size[0] * m] = 1.0;
    }
  } else {
    for (m = 0; m <= nv; m++) {
      d->data[m + d->size[0] * (m + absk)] = 1.0;
    }
  }
}

void f_diag(const emxArray_real_T *v, emxArray_real_T *d)
{
  int i;
  int loop_ub;
  int nv;
  nv = v->size[0];
  i = d->size[0] * d->size[1];
  d->size[0] = v->size[0];
  d->size[1] = v->size[0];
  emxEnsureCapacity_real_T(d, i);
  loop_ub = v->size[0] * v->size[0];
  for (i = 0; i < loop_ub; i++) {
    d->data[i] = 0.0;
  }

  for (loop_ub = 0; loop_ub < nv; loop_ub++) {
    d->data[loop_ub + d->size[0] * loop_ub] = 1.0;
  }
}

void g_diag(int K,const double v[K], double d[K*K])
{
  int j;
  memset(&d[0], 0, K*K * sizeof(double));
  for (j = 0; j < K; j++) {
    d[j + K * j] = v[j];
  }
}

/* End of code generation (diag.c) */
