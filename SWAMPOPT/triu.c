/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * triu.c
 *
 * Code generation for function 'triu'
 *
 */

/* Include files */
#include "triu.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void triu(emxArray_real_T *x)
{
  int i;
  int istart;
  int j;
  int jend;
  int m;
  m = x->size[0];
  if ((x->size[0] != 0) && (x->size[1] != 0) && (1 < x->size[0])) {
    istart = 2;
    if (x->size[0] - 2 < x->size[1] - 1) {
      jend = x->size[0] - 1;
    } else {
      jend = x->size[1];
    }

    for (j = 0; j < jend; j++) {
      for (i = istart; i <= m; i++) {
        x->data[(i + x->size[0] * j) - 1] = 0.0;
      }

      istart++;
    }
  }
}

/* End of code generation (triu.c) */
