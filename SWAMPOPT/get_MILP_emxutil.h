/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * get_MILP_emxutil.h
 *
 * Code generation for function 'get_MILP_emxutil'
 *
 */

#ifndef GET_MILP_EMXUTIL_H
#define GET_MILP_EMXUTIL_H

/* Include files */
#include "get_MILP_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  extern void emxEnsureCapacity_boolean_T(emxArray_boolean_T *emxArray, int
    oldNumel);
  extern void emxEnsureCapacity_char_T(emxArray_char_T *emxArray, int oldNumel);
  extern void emxEnsureCapacity_int32_T(emxArray_int32_T *emxArray, int oldNumel);
  extern void emxEnsureCapacity_int8_T(emxArray_int8_T *emxArray, int oldNumel);
  extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);
  extern void emxFreeMatrix_cell_wrap_6(cell_wrap_6 pMatrix[2]);
  extern void emxFreeStruct_cell_wrap_6(cell_wrap_6 *pStruct);
  extern void emxFreeStruct_sparse(d_sparse *pStruct);
  extern void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);
  extern void emxFree_char_T(emxArray_char_T **pEmxArray);
  extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
  extern void emxFree_int8_T(emxArray_int8_T **pEmxArray);
  extern void emxFree_real_T(emxArray_real_T **pEmxArray);
  extern void emxInitMatrix_cell_wrap_6(cell_wrap_6 pMatrix[2]);
  extern void emxInitStruct_cell_wrap_6(cell_wrap_6 *pStruct);
  extern void emxInitStruct_sparse(d_sparse *pStruct);
  extern void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int
    numDimensions);
  extern void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions);
  extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);
  extern void emxInit_int8_T(emxArray_int8_T **pEmxArray, int numDimensions);
  extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (get_MILP_emxutil.h) */
