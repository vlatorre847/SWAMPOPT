/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 */

/* Include files */
#include "repmat.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_repmat(const emxArray_real_T *a, double varargin_1, emxArray_real_T *b)
{
  int ibcol;
  int itilerow;
  int k;
  int nrows;
  int ntilerows;
  nrows = b->size[0];
  b->size[0] = a->size[0] * (int)varargin_1;
  emxEnsureCapacity_real_T(b, nrows);
  nrows = a->size[0];
  ntilerows = (int)varargin_1;
  for (itilerow = 0; itilerow < ntilerows; itilerow++) {
    ibcol = itilerow * nrows;
    for (k = 0; k < nrows; k++) {
      b->data[ibcol + k] = 1.0;
    }
  }
}

void c_repmat(const emxArray_real_T *a, double varargin_2, emxArray_real_T *b)
{
  int iacol_tmp;
  int ibmat;
  int ibtile;
  int jcol;
  int jtilecol;
  int k;
  int ncols;
  int nrows;
  int ntilecols;
  nrows = b->size[0] * b->size[1];
  b->size[0] = a->size[0];
  b->size[1] = a->size[1] * (int)varargin_2;
  emxEnsureCapacity_real_T(b, nrows);
  nrows = a->size[0];
  ncols = a->size[1];
  ntilecols = (int)varargin_2;
  for (jtilecol = 0; jtilecol < ntilecols; jtilecol++) {
    ibtile = jtilecol * (nrows * ncols) - 1;
    for (jcol = 0; jcol < ncols; jcol++) {
      iacol_tmp = jcol * nrows;
      ibmat = ibtile + iacol_tmp;
      for (k = 0; k < nrows; k++) {
        b->data[(ibmat + k) + 1] = a->data[iacol_tmp + k];
      }
    }
  }
}

void d_repmat(int l,const double a[l*l], double varargin_2, emxArray_real_T *b)
{
  int ibtile;
  int jcol;
  int jtilecol;
  int ntilecols;
  ntilecols = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = (l*l) * (int)varargin_2;
  emxEnsureCapacity_real_T(b, ntilecols);
  ntilecols = (int)varargin_2;
  for (jtilecol = 0; jtilecol < ntilecols; jtilecol++) {
    ibtile = jtilecol * (l*l);
    for (jcol = 0; jcol < (l*l); jcol++) {
      b->data[ibtile + jcol] = a[jcol];
    }
  }
}

void repmat(int l,const double a[l], double varargin_1, emxArray_real_T *b)
{
  int ibcol;
  int itilerow;
  int k;
  int ntilerows;
  ntilerows = b->size[0];
  b->size[0] = l * (int)varargin_1;
  emxEnsureCapacity_real_T(b, ntilerows);
  ntilerows = (int)varargin_1;
  for (itilerow = 0; itilerow < ntilerows; itilerow++) {
    ibcol = itilerow * l;
    for (k = 0; k < l; k++) {
      b->data[ibcol + k] = a[k];
    }
  }
}

/* End of code generation (repmat.c) */
