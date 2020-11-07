/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sparse1.c
 *
 * Code generation for function 'sparse1'
 *
 */

/* Include files */
#include "sparse1.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "locBsearch.h"
#include "parenAssign2D.h"
#include "rt_nonfinite.h"
#include <stddef.h>
#include <string.h>

/* Function Definitions */
void sparse_ctranspose(const emxArray_real_T *this_d, const emxArray_int32_T
  *this_colidx, const emxArray_int32_T *this_rowidx, int this_n, emxArray_real_T
  *y_d, emxArray_int32_T *y_colidx, emxArray_int32_T *y_rowidx, int *y_m)
{
  int c;
  int idx;
  int numalloc;
  int outridx;
  if (this_colidx->data[this_colidx->size[0] - 1] - 1 >= 1) {
    numalloc = this_colidx->data[this_colidx->size[0] - 1] - 2;
  } else {
    numalloc = 0;
  }

  idx = y_d->size[0];
  y_d->size[0] = numalloc + 1;
  emxEnsureCapacity_real_T(y_d, idx);
  for (idx = 0; idx <= numalloc; idx++) {
    y_d->data[idx] = 0.0;
  }

  idx = y_colidx->size[0];
  y_colidx->size[0] = 2;
  emxEnsureCapacity_int32_T(y_colidx, idx);
  idx = y_rowidx->size[0];
  y_rowidx->size[0] = numalloc + 1;
  emxEnsureCapacity_int32_T(y_rowidx, idx);
  for (idx = 0; idx <= numalloc; idx++) {
    y_rowidx->data[idx] = 0;
  }

  y_colidx->data[0] = 1;
  y_colidx->data[1] = 1;
  if (this_n != 0) {
    idx = y_colidx->size[0];
    y_colidx->size[0] = 2;
    emxEnsureCapacity_int32_T(y_colidx, idx);
    for (idx = 0; idx < 2; idx++) {
      y_colidx->data[idx] = 0;
    }

    idx = this_colidx->data[this_colidx->size[0] - 1];
    for (numalloc = 0; numalloc <= idx - 2; numalloc++) {
      y_colidx->data[this_rowidx->data[numalloc]]++;
    }

    y_colidx->data[0] = 1;
    y_colidx->data[1] += y_colidx->data[0];
    numalloc = -1;
    for (c = 0; c < this_n; c++) {
      for (idx = this_colidx->data[c]; idx < this_colidx->data[c + 1]; idx++) {
        outridx = numalloc + y_colidx->data[0];
        y_d->data[outridx] = this_d->data[idx - 1];
        y_rowidx->data[outridx] = c + 1;
        numalloc++;
      }
    }
  }

  *y_m = this_n;
}

void sparse_full(const emxArray_real_T *this_d, const emxArray_int32_T
                 *this_colidx, const emxArray_int32_T *this_rowidx, int this_m,
                 emxArray_real_T *y)
{
  int cend;
  int i;
  int idx;
  i = y->size[0] * y->size[1];
  y->size[0] = this_m;
  y->size[1] = 1;
  emxEnsureCapacity_real_T(y, i);
  for (i = 0; i < this_m; i++) {
    y->data[i] = 0.0;
  }

  cend = this_colidx->data[1] - 1;
  i = this_colidx->data[0];
  for (idx = i; idx <= cend; idx++) {
    y->data[this_rowidx->data[idx - 1] - 1] = this_d->data[idx - 1];
  }
}

void sparse_ne(const emxArray_real_T *a_d, const emxArray_int32_T *a_colidx,
               const emxArray_int32_T *a_rowidx, int a_m, int a_n,
               emxArray_boolean_T *s_d, emxArray_int32_T *s_colidx,
               emxArray_int32_T *s_rowidx, int *s_m, int *s_n)
{
  emxArray_boolean_T *tmpd;
  int c;
  int currRowIdx;
  int i;
  int numalloc;
  int ridx;
  boolean_T val;
  ridx = a_colidx->data[a_colidx->size[0] - 1];
  if (1 > a_colidx->data[a_colidx->size[0] - 1] - 1) {
    numalloc = 0;
  } else {
    numalloc = a_colidx->data[a_colidx->size[0] - 1] - 1;
  }

  emxInit_boolean_T(&tmpd, 1);
  i = tmpd->size[0];
  tmpd->size[0] = numalloc;
  emxEnsureCapacity_boolean_T(tmpd, i);
  for (i = 0; i < numalloc; i++) {
    tmpd->data[i] = (a_d->data[i] != 0.0);
  }

  if (a_colidx->data[a_colidx->size[0] - 1] - 1 >= 1) {
    numalloc = a_colidx->data[a_colidx->size[0] - 1] - 2;
  } else {
    numalloc = 0;
  }

  i = s_d->size[0];
  s_d->size[0] = numalloc + 1;
  emxEnsureCapacity_boolean_T(s_d, i);
  for (i = 0; i <= numalloc; i++) {
    s_d->data[i] = false;
  }

  i = s_rowidx->size[0];
  s_rowidx->size[0] = numalloc + 1;
  emxEnsureCapacity_int32_T(s_rowidx, i);
  for (i = 0; i <= numalloc; i++) {
    s_rowidx->data[i] = 0;
  }

  if (1 > a_colidx->data[a_colidx->size[0] - 1] - 1) {
    numalloc = 1;
  } else {
    numalloc = a_colidx->data[a_colidx->size[0] - 1];
  }

  for (i = 0; i <= numalloc - 2; i++) {
    s_rowidx->data[i] = a_rowidx->data[i];
  }

  i = s_colidx->size[0];
  s_colidx->size[0] = a_colidx->size[0];
  emxEnsureCapacity_int32_T(s_colidx, i);
  numalloc = a_colidx->size[0];
  for (i = 0; i < numalloc; i++) {
    s_colidx->data[i] = a_colidx->data[i];
  }

  for (numalloc = 0; numalloc <= ridx - 2; numalloc++) {
    s_d->data[numalloc] = tmpd->data[numalloc];
  }

  emxFree_boolean_T(&tmpd);
  numalloc = 1;
  i = a_colidx->size[0];
  for (c = 0; c <= i - 2; c++) {
    ridx = s_colidx->data[c];
    s_colidx->data[c] = numalloc;
    while (ridx < s_colidx->data[c + 1]) {
      currRowIdx = s_rowidx->data[ridx - 1];
      val = s_d->data[ridx - 1];
      ridx++;
      if (val) {
        s_d->data[numalloc - 1] = true;
        s_rowidx->data[numalloc - 1] = currRowIdx;
        numalloc++;
      }
    }
  }

  s_colidx->data[s_colidx->size[0] - 1] = numalloc;
  *s_m = a_m;
  *s_n = a_n;
}

void sparse_nonzeros(const emxArray_real_T *this_d, const emxArray_int32_T
                     *this_colidx, emxArray_real_T *y)
{
  int i;
  int loop_ub;
  if (1 > this_colidx->data[this_colidx->size[0] - 1] - 1) {
    loop_ub = 0;
  } else {
    loop_ub = this_colidx->data[this_colidx->size[0] - 1] - 1;
  }

  i = y->size[0];
  y->size[0] = loop_ub;
  emxEnsureCapacity_real_T(y, i);
  for (i = 0; i < loop_ub; i++) {
    y->data[i] = this_d->data[i];
  }
}

void sparse_parenAssign(d_sparse *this, int rhs_n, double varargin_1)
{
  double thisv;
  int cidx;
  int k;
  int nelem;
  int sn;
  int vidx;
  boolean_T found;
  if (rhs_n == 1) {
    sn = this->n;
    for (cidx = 0; cidx < sn; cidx++) {
      sparse_locBsearch(this->rowidx, (int)varargin_1, this->colidx->data[cidx],
                        this->colidx->data[cidx + 1], &vidx, &found);
      if (found) {
        thisv = this->d->data[vidx - 1];
      } else {
        thisv = 0.0;
      }

      if (!(thisv == 0.0)) {
        if (thisv == 0.0) {
          if (this->colidx->data[this->colidx->size[0] - 1] - 1 == this->maxnz)
          {
            b_realloc(this, this->colidx->data[this->colidx->size[0] - 1] + 9,
                      vidx, vidx + 1, this->colidx->data[this->colidx->size[0] -
                      1] - 1);
            this->rowidx->data[vidx] = (int)varargin_1;
            this->d->data[vidx] = 0.0;
          } else {
            nelem = (this->colidx->data[this->colidx->size[0] - 1] - vidx) - 1;
            if (nelem > 0) {
              memmove((void *)&this->rowidx->data[vidx + 1], (void *)
                      &this->rowidx->data[vidx], (unsigned int)((size_t)nelem *
                       sizeof(int)));
              memmove((void *)&this->d->data[vidx + 1], (void *)&this->d->
                      data[vidx], (unsigned int)((size_t)nelem * sizeof(double)));
            }

            this->d->data[vidx] = 0.0;
            this->rowidx->data[vidx] = (int)varargin_1;
          }

          vidx = cidx + 2;
          nelem = this->n + 1;
          for (k = vidx; k <= nelem; k++) {
            this->colidx->data[k - 1]++;
          }
        } else {
          nelem = (this->colidx->data[this->colidx->size[0] - 1] - vidx) - 1;
          if (nelem > 0) {
            memmove((void *)&this->rowidx->data[vidx - 1], (void *)&this->
                    rowidx->data[vidx], (unsigned int)((size_t)nelem * sizeof
                     (int)));
            memmove((void *)&this->d->data[vidx - 1], (void *)&this->d->
                    data[vidx], (unsigned int)((size_t)nelem * sizeof(double)));
          }

          vidx = cidx + 2;
          nelem = this->n + 1;
          for (k = vidx; k <= nelem; k++) {
            this->colidx->data[k - 1]--;
          }
        }
      }
    }
  } else {
    sn = this->n;
    for (cidx = 0; cidx < sn; cidx++) {
      sparse_locBsearch(this->rowidx, (int)varargin_1, this->colidx->data[cidx],
                        this->colidx->data[cidx + 1], &vidx, &found);
      if (found) {
        thisv = this->d->data[vidx - 1];
      } else {
        thisv = 0.0;
      }

      if (!(thisv == 0.0)) {
        if (thisv == 0.0) {
          if (this->colidx->data[this->colidx->size[0] - 1] - 1 == this->maxnz)
          {
            b_realloc(this, this->colidx->data[this->colidx->size[0] - 1] + 9,
                      vidx, vidx + 1, this->colidx->data[this->colidx->size[0] -
                      1] - 1);
            this->rowidx->data[vidx] = (int)varargin_1;
            this->d->data[vidx] = 0.0;
          } else {
            nelem = (this->colidx->data[this->colidx->size[0] - 1] - vidx) - 1;
            if (nelem > 0) {
              memmove((void *)&this->rowidx->data[vidx + 1], (void *)
                      &this->rowidx->data[vidx], (unsigned int)((size_t)nelem *
                       sizeof(int)));
              memmove((void *)&this->d->data[vidx + 1], (void *)&this->d->
                      data[vidx], (unsigned int)((size_t)nelem * sizeof(double)));
            }

            this->d->data[vidx] = 0.0;
            this->rowidx->data[vidx] = (int)varargin_1;
          }

          vidx = cidx + 2;
          nelem = this->n + 1;
          for (k = vidx; k <= nelem; k++) {
            this->colidx->data[k - 1]++;
          }
        } else {
          nelem = (this->colidx->data[this->colidx->size[0] - 1] - vidx) - 1;
          if (nelem > 0) {
            memmove((void *)&this->rowidx->data[vidx - 1], (void *)&this->
                    rowidx->data[vidx], (unsigned int)((size_t)nelem * sizeof
                     (int)));
            memmove((void *)&this->d->data[vidx - 1], (void *)&this->d->
                    data[vidx], (unsigned int)((size_t)nelem * sizeof(double)));
          }

          vidx = cidx + 2;
          nelem = this->n + 1;
          for (k = vidx; k <= nelem; k++) {
            this->colidx->data[k - 1]--;
          }
        }
      }
    }
  }
}

void sparse_parenReference(const emxArray_real_T *this_d, const emxArray_int32_T
  *this_colidx, const emxArray_int32_T *this_rowidx, int this_n, const
  emxArray_real_T *varargin_1, emxArray_real_T *s_d, emxArray_int32_T *s_colidx,
  emxArray_int32_T *s_rowidx, int *s_m, int *s_n, int *s_maxnz)
{
  double s_d_tmp;
  int cidx;
  int colNnz;
  int i;
  int i1;
  int idx;
  int k;
  int ridx;
  int sm;
  boolean_T found;
  sm = varargin_1->size[1];
  s_d->size[0] = 0;
  s_rowidx->size[0] = 0;
  i = s_colidx->size[0];
  s_colidx->size[0] = this_n + 1;
  emxEnsureCapacity_int32_T(s_colidx, i);
  for (i = 0; i <= this_n; i++) {
    s_colidx->data[i] = 0;
  }

  s_colidx->data[0] = 1;
  colNnz = 1;
  k = 0;
  for (cidx = 0; cidx < this_n; cidx++) {
    for (ridx = 0; ridx < sm; ridx++) {
      sparse_locBsearch(this_rowidx, (int)varargin_1->data[ridx],
                        this_colidx->data[cidx], this_colidx->data[cidx + 1],
                        &idx, &found);
      if (found) {
        i = s_d->size[0];
        i1 = s_d->size[0];
        s_d->size[0]++;
        emxEnsureCapacity_real_T(s_d, i1);
        s_d_tmp = this_d->data[idx - 1];
        s_d->data[i] = s_d_tmp;
        i = s_rowidx->size[0];
        i1 = s_rowidx->size[0];
        s_rowidx->size[0]++;
        emxEnsureCapacity_int32_T(s_rowidx, i1);
        s_rowidx->data[i] = ridx + 1;
        s_d->data[k] = s_d_tmp;
        s_rowidx->data[k] = ridx + 1;
        k++;
        colNnz++;
      }
    }

    s_colidx->data[cidx + 1] = colNnz;
  }

  if (s_colidx->data[s_colidx->size[0] - 1] - 1 == 0) {
    i = s_rowidx->size[0];
    s_rowidx->size[0] = 1;
    emxEnsureCapacity_int32_T(s_rowidx, i);
    s_rowidx->data[0] = 1;
    i = s_d->size[0];
    s_d->size[0] = 1;
    emxEnsureCapacity_real_T(s_d, i);
    s_d->data[0] = 0.0;
  }

  *s_m = varargin_1->size[1];
  if (s_colidx->data[s_colidx->size[0] - 1] - 1 >= 1) {
    *s_maxnz = s_colidx->data[s_colidx->size[0] - 1] - 1;
  } else {
    *s_maxnz = 1;
  }

  *s_n = this_n;
}

/* End of code generation (sparse1.c) */
