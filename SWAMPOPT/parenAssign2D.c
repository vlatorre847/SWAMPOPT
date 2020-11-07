/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * parenAssign2D.c
 *
 * Code generation for function 'parenAssign2D'
 *
 */

/* Include files */
#include "parenAssign2D.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_realloc(d_sparse *this, int numAllocRequested, int ub1, int lb2, int ub2)
{
  emxArray_int32_T *rowidxt;
  emxArray_real_T *dt;
  int highOrderA;
  int highOrderB;
  int lowOrderB;
  int overflow;
  int partialResults_idx_0_tmp;
  int partialResults_idx_1;
  int tmp;
  emxInit_int32_T(&rowidxt, 1);
  highOrderA = rowidxt->size[0];
  rowidxt->size[0] = this->rowidx->size[0];
  emxEnsureCapacity_int32_T(rowidxt, highOrderA);
  partialResults_idx_1 = this->rowidx->size[0];
  for (highOrderA = 0; highOrderA < partialResults_idx_1; highOrderA++) {
    rowidxt->data[highOrderA] = this->rowidx->data[highOrderA];
  }

  emxInit_real_T(&dt, 1);
  highOrderA = dt->size[0];
  dt->size[0] = this->d->size[0];
  emxEnsureCapacity_real_T(dt, highOrderA);
  partialResults_idx_1 = this->d->size[0];
  for (highOrderA = 0; highOrderA < partialResults_idx_1; highOrderA++) {
    dt->data[highOrderA] = this->d->data[highOrderA];
  }

  highOrderA = this->m >> 16;
  partialResults_idx_1 = this->m & 65535;
  highOrderB = this->n >> 16;
  lowOrderB = this->n & 65535;
  partialResults_idx_0_tmp = partialResults_idx_1 * lowOrderB;
  tmp = partialResults_idx_1 * highOrderB;
  partialResults_idx_1 = tmp << 16;
  overflow = tmp >> 16;
  if (overflow <= 0) {
    tmp = highOrderA * lowOrderB;
    overflow += tmp >> 16;
    if (overflow <= 0) {
      overflow += highOrderA * highOrderB;
      if (overflow <= 0) {
        if (partialResults_idx_0_tmp > MAX_int32_T - partialResults_idx_1) {
          partialResults_idx_1 = (partialResults_idx_0_tmp +
            partialResults_idx_1) - MAX_int32_T;
          overflow++;
        } else {
          partialResults_idx_1 += partialResults_idx_0_tmp;
        }

        if (partialResults_idx_1 > MAX_int32_T - (tmp << 16)) {
          overflow++;
        }
      }
    }
  }

  if (overflow == 0) {
    partialResults_idx_1 = this->m * this->n;
    if (numAllocRequested <= partialResults_idx_1) {
      partialResults_idx_1 = numAllocRequested;
    }

    if (1 >= partialResults_idx_1) {
      partialResults_idx_1 = 1;
    }
  } else if (1 >= numAllocRequested) {
    partialResults_idx_1 = 1;
  } else {
    partialResults_idx_1 = numAllocRequested;
  }

  highOrderA = this->rowidx->size[0];
  this->rowidx->size[0] = partialResults_idx_1;
  emxEnsureCapacity_int32_T(this->rowidx, highOrderA);
  for (highOrderA = 0; highOrderA < partialResults_idx_1; highOrderA++) {
    this->rowidx->data[highOrderA] = 0;
  }

  highOrderA = this->d->size[0];
  this->d->size[0] = partialResults_idx_1;
  emxEnsureCapacity_real_T(this->d, highOrderA);
  for (highOrderA = 0; highOrderA < partialResults_idx_1; highOrderA++) {
    this->d->data[highOrderA] = 0.0;
  }

  this->maxnz = partialResults_idx_1;
  for (partialResults_idx_1 = 0; partialResults_idx_1 < ub1;
       partialResults_idx_1++) {
    this->rowidx->data[partialResults_idx_1] = rowidxt->
      data[partialResults_idx_1];
    this->d->data[partialResults_idx_1] = dt->data[partialResults_idx_1];
  }

  for (partialResults_idx_1 = lb2; partialResults_idx_1 <= ub2;
       partialResults_idx_1++) {
    this->rowidx->data[partialResults_idx_1] = rowidxt->
      data[partialResults_idx_1 - 1];
    this->d->data[partialResults_idx_1] = dt->data[partialResults_idx_1 - 1];
  }

  emxFree_real_T(&dt);
  emxFree_int32_T(&rowidxt);
}

/* End of code generation (parenAssign2D.c) */
