/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * get_MILP.h
 *
 * Code generation for function 'get_MILP'
 *
 */

#ifndef GET_MILP_H
#define GET_MILP_H

/* Include files */
#include "get_MILP_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  extern void get_MILP(int N, int M, int K, int l, int tp_size,
                       int tpi_size, int total_delay, const double V0[total_delay], const double R0[l],
                       const double H0[l], const int tp[tp_size], const int
                       tpi[tpi_size], const double c[l], const double r[N], const
                       double b_gamma[l], const double rho[l], const double psi
                       [l*l], const double delay[l], double eps_min, const int
                       Ki[K], const int Ii[l], const double q[K], const
                       int s[K], const int d[K], const double eps[K],
                       const double alfa[K], const double beta[K], const
                       double Dt[K], const double Dv[K], const double j_funct
                       [4], emxArray_real_T *f, emxArray_real_T *A_beg,
                       emxArray_real_T *A_nnz, emxArray_int32_T *A_indx,
                       emxArray_real_T *A_val, emxArray_char_T *A_sense,
                       emxArray_real_T *b, emxArray_real_T *lb, emxArray_real_T *ub,
                       emxArray_char_T *ctype, emxArray_real_T *A_j, double *A_r);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (get_MILP.h) */
