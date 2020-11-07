/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * repmat1.c
 *
 * Code generation for function 'repmat1'
 *
 */

/* Include files */
#include "repmat1.h"
#include "fillIn.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void sparse_repmat(const emxArray_real_T *A_d, const emxArray_int32_T *A_colidx,
                   const emxArray_int32_T *A_rowidx, int A_m, int A_n, double
                   varargin_2, d_sparse *B)
{
  int b_i;
  int ci;
  int i;
  int i1;
  int i2;
  int i3;
  int numalloc;
  int numelThisCol;
  numelThisCol = A_n * (int)varargin_2;
  numalloc = (A_colidx->data[A_colidx->size[0] - 1] - 1) * (int)varargin_2;
  B->m = A_m;
  B->n = numelThisCol;
  if (numalloc < 1) {
    numalloc = 1;
  }

  i = B->d->size[0];
  B->d->size[0] = numalloc;
  emxEnsureCapacity_real_T(B->d, i);
  for (i = 0; i < numalloc; i++) {
    B->d->data[i] = 0.0;
  }

  B->maxnz = numalloc;
  i = B->colidx->size[0];
  B->colidx->size[0] = numelThisCol + 1;
  emxEnsureCapacity_int32_T(B->colidx, i);
  B->colidx->data[0] = 1;
  i = B->rowidx->size[0];
  B->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(B->rowidx, i);
  for (i = 0; i < numalloc; i++) {
    B->rowidx->data[i] = 0;
  }

  for (numalloc = 0; numalloc < numelThisCol; numalloc++) {
    B->colidx->data[numalloc + 1] = 1;
  }

  sparse_fillIn(B);
  if ((B->m != 0) && (B->n != 0)) {
    B->colidx->data[0] = 1;
    numalloc = 0;
    for (ci = 0; ci < A_n; ci++) {
      i = A_colidx->data[ci];
      i1 = A_colidx->data[ci + 1];
      B->colidx->data[ci + 1] = (B->colidx->data[ci] + i1) - A_colidx->data[ci];
      numelThisCol = i1 - A_colidx->data[ci];
      if (numelThisCol > 0) {
        i1 = numelThisCol - 1;
        for (b_i = 0; b_i <= i1; b_i++) {
          i2 = (i + b_i) - 1;
          i3 = b_i + numalloc;
          B->d->data[i3] = A_d->data[i2];
          B->rowidx->data[i3] = A_rowidx->data[i2];
        }

        numalloc += numelThisCol;
      }
    }

    numalloc = B->colidx->data[A_n] - 1;
    i = (int)varargin_2;
    for (numelThisCol = 0; numelThisCol <= i - 2; numelThisCol++) {
      for (b_i = 0; b_i < numalloc; b_i++) {
        i1 = b_i + numalloc * (numelThisCol + 1);
        B->rowidx->data[i1] = B->rowidx->data[b_i];
        B->d->data[i1] = B->d->data[b_i];
      }

      for (b_i = 0; b_i < A_n; b_i++) {
        B->colidx->data[A_n * (numelThisCol + 1) + b_i] = B->colidx->data[b_i] +
          (numelThisCol + 1) * numalloc;
      }
    }

    B->colidx->data[B->colidx->size[0] - 1] = numalloc * (int)varargin_2 + 1;
  }
}

/* End of code generation (repmat1.c) */
