/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * heapsort.c
 *
 * Code generation for function 'heapsort'
 *
 */

/* Include files */
#include "heapsort.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void heapify(emxArray_int32_T *x, int idx, int xstart, int xend, const
                    cell_wrap_6 cmp_tunableEnvironment[2]);

/* Function Definitions */
static void heapify(emxArray_int32_T *x, int idx, int xstart, int xend, const
                    cell_wrap_6 cmp_tunableEnvironment[2])
{
  int cmpIdx;
  int extremum;
  int extremumIdx;
  int i;
  int i1;
  int i2;
  int leftIdx;
  int xcmp;
  boolean_T changed;
  boolean_T exitg1;
  boolean_T varargout_1;
  changed = true;
  extremumIdx = (idx + xstart) - 2;
  leftIdx = ((idx << 1) + xstart) - 2;
  exitg1 = false;
  while ((!exitg1) && (leftIdx + 1 < xend)) {
    changed = false;
    extremum = x->data[extremumIdx];
    cmpIdx = leftIdx;
    xcmp = x->data[leftIdx];
    i = x->data[leftIdx + 1] - 1;
    i1 = cmp_tunableEnvironment[0].f1->data[x->data[leftIdx] - 1];
    i2 = cmp_tunableEnvironment[0].f1->data[i];
    if (i1 < i2) {
      varargout_1 = true;
    } else if (i1 == i2) {
      varargout_1 = (cmp_tunableEnvironment[1].f1->data[x->data[leftIdx] - 1] <
                     cmp_tunableEnvironment[1].f1->data[i]);
    } else {
      varargout_1 = false;
    }

    if (varargout_1) {
      cmpIdx = leftIdx + 1;
      xcmp = x->data[leftIdx + 1];
    }

    i = cmp_tunableEnvironment[0].f1->data[x->data[extremumIdx] - 1];
    i1 = cmp_tunableEnvironment[0].f1->data[xcmp - 1];
    if (i < i1) {
      varargout_1 = true;
    } else if (i == i1) {
      varargout_1 = (cmp_tunableEnvironment[1].f1->data[x->data[extremumIdx] - 1]
                     < cmp_tunableEnvironment[1].f1->data[xcmp - 1]);
    } else {
      varargout_1 = false;
    }

    if (varargout_1) {
      x->data[extremumIdx] = xcmp;
      x->data[cmpIdx] = extremum;
      extremumIdx = cmpIdx;
      leftIdx = ((((cmpIdx - xstart) + 2) << 1) + xstart) - 2;
      changed = true;
    } else {
      exitg1 = true;
    }
  }

  if (changed && (leftIdx + 1 <= xend)) {
    extremum = x->data[extremumIdx];
    i = cmp_tunableEnvironment[0].f1->data[x->data[extremumIdx] - 1];
    i1 = cmp_tunableEnvironment[0].f1->data[x->data[leftIdx] - 1];
    if (i < i1) {
      varargout_1 = true;
    } else if (i == i1) {
      varargout_1 = (cmp_tunableEnvironment[1].f1->data[x->data[extremumIdx] - 1]
                     < cmp_tunableEnvironment[1].f1->data[x->data[leftIdx] - 1]);
    } else {
      varargout_1 = false;
    }

    if (varargout_1) {
      x->data[extremumIdx] = x->data[leftIdx];
      x->data[leftIdx] = extremum;
    }
  }
}

void b_heapsort(emxArray_int32_T *x, int xstart, int xend, const cell_wrap_6
                cmp_tunableEnvironment[2])
{
  int k;
  int n;
  int t;
  n = (xend - xstart) - 1;
  for (t = n + 2; t >= 1; t--) {
    heapify(x, t, xstart, xend, cmp_tunableEnvironment);
  }

  for (k = 0; k <= n; k++) {
    t = x->data[xend - 1];
    x->data[xend - 1] = x->data[xstart - 1];
    x->data[xstart - 1] = t;
    xend--;
    heapify(x, 1, xstart, xend, cmp_tunableEnvironment);
  }
}

/* End of code generation (heapsort.c) */
