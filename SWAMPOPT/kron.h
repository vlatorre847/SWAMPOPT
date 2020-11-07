/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * kron.h
 *
 * Code generation for function 'kron'
 *
 */

#ifndef KRON_H
#define KRON_H

/* Include files */
#include "get_MILP_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  void b_kron(int l,const emxArray_real_T *A, const double B[l*l], emxArray_real_T *K);
  void c_kron(const emxArray_real_T *A, const emxArray_real_T *B,
              emxArray_real_T *K);
  void d_kron(int l,const emxArray_real_T *A, const double B[l*l], emxArray_real_T *K);
  void e_kron(const emxArray_real_T *A, const emxArray_real_T *B,
              emxArray_real_T *K);
  void f_kron(int KK, const emxArray_real_T *A, const double B[KK*KK], emxArray_real_T *K);
  void g_kron(int l,const emxArray_real_T *A, const double B[l*l], emxArray_real_T *K);
  void kron(const emxArray_real_T *A, const emxArray_real_T *B, emxArray_real_T *
            K);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (kron.h) */
