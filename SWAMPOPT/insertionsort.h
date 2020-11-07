/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * insertionsort.h
 *
 * Code generation for function 'insertionsort'
 *
 */

#ifndef INSERTIONSORT_H
#define INSERTIONSORT_H

/* Include files */
#include "get_MILP_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  void insertionsort(emxArray_int32_T *x, int xstart, int xend, const
                     cell_wrap_6 cmp_tunableEnvironment[2]);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (insertionsort.h) */
