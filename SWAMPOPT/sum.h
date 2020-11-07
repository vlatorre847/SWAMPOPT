/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum.h
 *
 * Code generation for function 'sum'
 *
 */

#ifndef SUM_H
#define SUM_H

/* Include files */
#include "get_MILP_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  void b_sum(const emxArray_boolean_T *x_d, const emxArray_int32_T *x_colidx,
             int x_m, int x_n, emxArray_real_T *y_d, emxArray_int32_T *y_colidx,
             emxArray_int32_T *y_rowidx, int *y_n);
  double sum(const double x_data[], const int x_size[1]);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (sum.h) */
