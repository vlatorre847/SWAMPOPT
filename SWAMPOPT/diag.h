/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * diag.h
 *
 * Code generation for function 'diag'
 *
 */

#ifndef DIAG_H
#define DIAG_H

/* Include files */
#include "get_MILP_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  void b_diag(int l,const double v[l], double d[l*l]);
  void c_diag(const emxArray_real_T *v, emxArray_real_T *d);
  void d_diag(const emxArray_real_T *v, emxArray_real_T *d);
  void diag(const emxArray_real_T *v, emxArray_real_T *d);
  void e_diag(const emxArray_real_T *v, double K, emxArray_real_T *d);
  void f_diag(const emxArray_real_T *v, emxArray_real_T *d);
  void g_diag(int K,const double v[K], double d[K*K]);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (diag.h) */
