/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * insertionsort.c
 *
 * Code generation for function 'insertionsort'
 *
 */

/* Include files */
#include "insertionsort.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void insertionsort(emxArray_int32_T *x, int xstart, int xend, const cell_wrap_6
                   cmp_tunableEnvironment[2])
{
  int i;
  int idx;
  int k;
  int xc;
  boolean_T exitg1;
  boolean_T varargout_1;
  i = xstart + 1;
  for (k = i; k <= xend; k++) {
    xc = x->data[k - 1] - 1;
    idx = k - 2;
    exitg1 = false;
    while ((!exitg1) && (idx + 1 >= xstart)) {
      if (cmp_tunableEnvironment[0].f1->data[xc] < cmp_tunableEnvironment[0].
          f1->data[x->data[idx] - 1]) {
        varargout_1 = true;
      } else if (cmp_tunableEnvironment[0].f1->data[xc] ==
                 cmp_tunableEnvironment[0].f1->data[x->data[idx] - 1]) {
        varargout_1 = (cmp_tunableEnvironment[1].f1->data[xc] <
                       cmp_tunableEnvironment[1].f1->data[x->data[idx] - 1]);
      } else {
        varargout_1 = false;
      }

      if (varargout_1) {
        x->data[idx + 1] = x->data[idx];
        idx--;
      } else {
        exitg1 = true;
      }
    }

    x->data[idx + 1] = xc + 1;
  }
}

/* End of code generation (insertionsort.c) */
