/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * repelem.c
 *
 * Code generation for function 'repelem'
 *
 */

/* Include files */
#include "repelem.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_repelem(double varargin_1, emxArray_char_T *y)
{
  int i;
  int idx;
  int j;
  i = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)varargin_1;
  emxEnsureCapacity_char_T(y, i);
  idx = -1;
  i = (int)varargin_1;
  for (j = 0; j < i; j++) {
    idx++;
    y->data[idx] = 'E';
  }
}

void repelem(double varargin_1, emxArray_char_T *y)
{
  int i;
  int idx;
  int j;
  i = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)varargin_1;
  emxEnsureCapacity_char_T(y, i);
  idx = -1;
  i = (int)varargin_1;
  for (j = 0; j < i; j++) {
    idx++;
    y->data[idx] = 'L';
  }
}

/* End of code generation (repelem.c) */
