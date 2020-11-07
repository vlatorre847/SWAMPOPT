/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eye.c
 *
 * Code generation for function 'eye'
 *
 */

/* Include files */
#include "eye.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void eye(double varargin_1, emxArray_real_T *b_I)
{
  double t;
  int i;
  int loop_ub;
  int m;
  if (varargin_1 < 0.0) {
    t = 0.0;
  } else {
    t = varargin_1;
  }

  m = (int)t;
  i = b_I->size[0] * b_I->size[1];
  b_I->size[0] = (int)t;
  b_I->size[1] = (int)t;
  emxEnsureCapacity_real_T(b_I, i);
  loop_ub = (int)t * (int)t;
  for (i = 0; i < loop_ub; i++) {
    b_I->data[i] = 0.0;
  }

  if ((int)t > 0) {
    for (loop_ub = 0; loop_ub < m; loop_ub++) {
      b_I->data[loop_ub + b_I->size[0] * loop_ub] = 1.0;
    }
  }
}

/* End of code generation (eye.c) */
