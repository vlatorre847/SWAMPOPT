/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * kron.c
 *
 * Code generation for function 'kron'
 *
 */

/* Include files */
#include "kron.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_kron(int l,const emxArray_real_T *A, const double B[l], emxArray_real_T *K)
{
  int b_j1;
  int i1;
  int i2;
  int j2;
  int kidx;
  int ma;
  int na;
  ma = A->size[0];
  na = A->size[1];
  kidx = K->size[0] * K->size[1];
  K->size[0] = A->size[0] * l;
  K->size[1] = A->size[1] * l;
  emxEnsureCapacity_real_T(K, kidx);
  kidx = -1;
  for (b_j1 = 0; b_j1 < na; b_j1++) {
    for (j2 = 0; j2 < l; j2++) {
      for (i1 = 0; i1 < ma; i1++) {
        for (i2 = 0; i2 < l; i2++) {
          kidx++;
          K->data[kidx] = A->data[i1 + A->size[0] * b_j1] * B[i2 + l * j2];
        }
      }
    }
  }
}

void c_kron(const emxArray_real_T *A, const emxArray_real_T *B, emxArray_real_T *
            K)
{
  int b_j1;
  int i1;
  int i2;
  int j2;
  int kidx;
  int ma;
  int mb;
  int na;
  int nb;
  ma = A->size[0];
  na = A->size[1];
  mb = B->size[0];
  nb = B->size[1];
  kidx = K->size[0] * K->size[1];
  K->size[0] = A->size[0] * B->size[0];
  K->size[1] = A->size[1] * B->size[1];
  emxEnsureCapacity_real_T(K, kidx);
  kidx = -1;
  for (b_j1 = 0; b_j1 < na; b_j1++) {
    for (j2 = 0; j2 < nb; j2++) {
      for (i1 = 0; i1 < ma; i1++) {
        for (i2 = 0; i2 < mb; i2++) {
          kidx++;
          K->data[kidx] = A->data[i1 + A->size[0] * b_j1] * B->data[i2 + B->
            size[0] * j2];
        }
      }
    }
  }
}

void d_kron(int l,const emxArray_real_T *A, const double B[l*l], emxArray_real_T *K)
{
  int b_j1;
  int i1;
  int j2;
  int kidx;
  int ma;
  int na;
  ma = A->size[0];
  na = A->size[1];
  kidx = K->size[0] * K->size[1];
  K->size[0] = A->size[0];
  K->size[1] = A->size[1] * l*l;
  emxEnsureCapacity_real_T(K, kidx);
  kidx = -1;
  for (b_j1 = 0; b_j1 < na; b_j1++) {
    for (j2 = 0; j2 < l*l; j2++) {
      for (i1 = 0; i1 < ma; i1++) {
        kidx++;
        K->data[kidx] = A->data[i1 + A->size[0] * b_j1] * B[j2];
      }
    }
  }
}

void e_kron(const emxArray_real_T *A, const emxArray_real_T *B, emxArray_real_T *
            K)
{
  int b_j1;
  int i2;
  int j2;
  int kidx;
  int mb;
  int na;
  int nb;
  na = A->size[1];
  mb = B->size[0];
  nb = B->size[1];
  kidx = K->size[0] * K->size[1];
  K->size[0] = B->size[0];
  K->size[1] = A->size[1] * B->size[1];
  emxEnsureCapacity_real_T(K, kidx);
  kidx = -1;
  for (b_j1 = 0; b_j1 < na; b_j1++) {
    for (j2 = 0; j2 < nb; j2++) {
      for (i2 = 0; i2 < mb; i2++) {
        kidx++;
        K->data[kidx] = A->data[b_j1] * B->data[i2 + B->size[0] * j2];
      }
    }
  }
}

void f_kron(int KK,const emxArray_real_T *A, const double B[KK*KK], emxArray_real_T *K)
{
  int b_j1;
  int i2;
  int j2;
  int kidx;
  int na;
  na = A->size[1];
  kidx = K->size[0] * K->size[1];
  K->size[0] = KK;
  K->size[1] = A->size[1] * KK;
  emxEnsureCapacity_real_T(K, kidx);
  kidx = -1;
  for (b_j1 = 0; b_j1 < na; b_j1++) {
    for (j2 = 0; j2 < KK; j2++) {
      for (i2 = 0; i2 < KK; i2++) {
        kidx++;
        K->data[kidx] = B[i2 + KK * j2];
      }
    }
  }
}

void g_kron(int l,const emxArray_real_T *A, const double B[l*l], emxArray_real_T *K)
{
  int i1;
  int i2;
  int kidx;
  int ma;
  ma = A->size[0];
  kidx = K->size[0];
  K->size[0] = A->size[0] * l*l;
  emxEnsureCapacity_real_T(K, kidx);
  kidx = -1;
  for (i1 = 0; i1 < ma; i1++) {
    for (i2 = 0; i2 < l*l; i2++) {
      kidx++;
      K->data[kidx] = B[i2];
    }
  }
}

void kron(const emxArray_real_T *A, const emxArray_real_T *B, emxArray_real_T *K)
{
  int b_j1;
  int i1;
  int j2;
  int kidx;
  int ma;
  int na;
  int nb;
  ma = A->size[0];
  na = A->size[1];
  nb = B->size[1];
  kidx = K->size[0] * K->size[1];
  K->size[0] = A->size[0];
  K->size[1] = A->size[1] * B->size[1];
  emxEnsureCapacity_real_T(K, kidx);
  kidx = -1;
  for (b_j1 = 0; b_j1 < na; b_j1++) {
    for (j2 = 0; j2 < nb; j2++) {
      for (i1 = 0; i1 < ma; i1++) {
        kidx++;
        K->data[kidx] = A->data[i1 + A->size[0] * b_j1] * B->data[j2];
      }
    }
  }
}

/* End of code generation (kron.c) */
