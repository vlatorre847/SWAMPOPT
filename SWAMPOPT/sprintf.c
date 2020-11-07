/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sprintf.c
 *
 * Code generation for function 'sprintf'
 *
 */

/* Include files */
#include "sprintf.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_types.h"
#include "rt_nonfinite.h"
#include <stddef.h>
#include <stdio.h>

/* Function Definitions */
void b_sprintf(const emxArray_char_T *varargin_1, const emxArray_char_T
               *varargin_2, const emxArray_char_T *varargin_3, const
               emxArray_char_T *varargin_4, const emxArray_char_T *varargin_5,
               const emxArray_char_T *varargin_6, const emxArray_char_T
               *varargin_7, const emxArray_char_T *varargin_8, const
               emxArray_char_T *varargin_9, const emxArray_char_T *varargin_10,
               emxArray_char_T *str)
{
  emxArray_char_T *b_varargin_1;
  emxArray_char_T *b_varargin_10;
  emxArray_char_T *b_varargin_2;
  emxArray_char_T *b_varargin_3;
  emxArray_char_T *b_varargin_4;
  emxArray_char_T *b_varargin_5;
  emxArray_char_T *b_varargin_6;
  emxArray_char_T *b_varargin_7;
  emxArray_char_T *b_varargin_8;
  emxArray_char_T *b_varargin_9;
  emxArray_char_T *c_varargin_1;
  emxArray_char_T *c_varargin_10;
  emxArray_char_T *c_varargin_2;
  emxArray_char_T *c_varargin_3;
  emxArray_char_T *c_varargin_4;
  emxArray_char_T *c_varargin_5;
  emxArray_char_T *c_varargin_6;
  emxArray_char_T *c_varargin_7;
  emxArray_char_T *c_varargin_8;
  emxArray_char_T *c_varargin_9;
  int i;
  int nbytes;
  emxInit_char_T(&b_varargin_1, 2);
  i = b_varargin_1->size[0] * b_varargin_1->size[1];
  b_varargin_1->size[0] = 1;
  b_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_1, i);
  nbytes = varargin_1->size[1];
  for (i = 0; i < nbytes; i++) {
    b_varargin_1->data[i] = varargin_1->data[i];
  }

  emxInit_char_T(&b_varargin_2, 2);
  b_varargin_1->data[varargin_1->size[1]] = '\x00';
  i = b_varargin_2->size[0] * b_varargin_2->size[1];
  b_varargin_2->size[0] = 1;
  b_varargin_2->size[1] = varargin_2->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_2, i);
  nbytes = varargin_2->size[1];
  for (i = 0; i < nbytes; i++) {
    b_varargin_2->data[i] = varargin_2->data[i];
  }

  emxInit_char_T(&b_varargin_3, 2);
  b_varargin_2->data[varargin_2->size[1]] = '\x00';
  i = b_varargin_3->size[0] * b_varargin_3->size[1];
  b_varargin_3->size[0] = 1;
  b_varargin_3->size[1] = varargin_3->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_3, i);
  nbytes = varargin_3->size[1];
  for (i = 0; i < nbytes; i++) {
    b_varargin_3->data[i] = varargin_3->data[i];
  }

  emxInit_char_T(&b_varargin_4, 2);
  b_varargin_3->data[varargin_3->size[1]] = '\x00';
  i = b_varargin_4->size[0] * b_varargin_4->size[1];
  b_varargin_4->size[0] = 1;
  b_varargin_4->size[1] = varargin_4->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_4, i);
  nbytes = varargin_4->size[1];
  for (i = 0; i < nbytes; i++) {
    b_varargin_4->data[i] = varargin_4->data[i];
  }

  emxInit_char_T(&b_varargin_5, 2);
  b_varargin_4->data[varargin_4->size[1]] = '\x00';
  i = b_varargin_5->size[0] * b_varargin_5->size[1];
  b_varargin_5->size[0] = 1;
  b_varargin_5->size[1] = varargin_5->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_5, i);
  nbytes = varargin_5->size[1];
  for (i = 0; i < nbytes; i++) {
    b_varargin_5->data[i] = varargin_5->data[i];
  }

  emxInit_char_T(&b_varargin_6, 2);
  b_varargin_5->data[varargin_5->size[1]] = '\x00';
  i = b_varargin_6->size[0] * b_varargin_6->size[1];
  b_varargin_6->size[0] = 1;
  b_varargin_6->size[1] = varargin_6->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_6, i);
  nbytes = varargin_6->size[1];
  for (i = 0; i < nbytes; i++) {
    b_varargin_6->data[i] = varargin_6->data[i];
  }

  emxInit_char_T(&b_varargin_7, 2);
  b_varargin_6->data[varargin_6->size[1]] = '\x00';
  i = b_varargin_7->size[0] * b_varargin_7->size[1];
  b_varargin_7->size[0] = 1;
  b_varargin_7->size[1] = varargin_7->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_7, i);
  nbytes = varargin_7->size[1];
  for (i = 0; i < nbytes; i++) {
    b_varargin_7->data[i] = varargin_7->data[i];
  }

  emxInit_char_T(&b_varargin_8, 2);
  b_varargin_7->data[varargin_7->size[1]] = '\x00';
  i = b_varargin_8->size[0] * b_varargin_8->size[1];
  b_varargin_8->size[0] = 1;
  b_varargin_8->size[1] = varargin_8->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_8, i);
  nbytes = varargin_8->size[1];
  for (i = 0; i < nbytes; i++) {
    b_varargin_8->data[i] = varargin_8->data[i];
  }

  emxInit_char_T(&b_varargin_9, 2);
  b_varargin_8->data[varargin_8->size[1]] = '\x00';
  i = b_varargin_9->size[0] * b_varargin_9->size[1];
  b_varargin_9->size[0] = 1;
  b_varargin_9->size[1] = varargin_9->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_9, i);
  nbytes = varargin_9->size[1];
  for (i = 0; i < nbytes; i++) {
    b_varargin_9->data[i] = varargin_9->data[i];
  }

  emxInit_char_T(&b_varargin_10, 2);
  b_varargin_9->data[varargin_9->size[1]] = '\x00';
  i = b_varargin_10->size[0] * b_varargin_10->size[1];
  b_varargin_10->size[0] = 1;
  b_varargin_10->size[1] = varargin_10->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_10, i);
  nbytes = varargin_10->size[1];
  for (i = 0; i < nbytes; i++) {
    b_varargin_10->data[i] = varargin_10->data[i];
  }

  emxInit_char_T(&c_varargin_1, 2);
  b_varargin_10->data[varargin_10->size[1]] = '\x00';
  i = c_varargin_1->size[0] * c_varargin_1->size[1];
  c_varargin_1->size[0] = 1;
  c_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_1, i);
  nbytes = varargin_1->size[1];
  for (i = 0; i < nbytes; i++) {
    c_varargin_1->data[i] = varargin_1->data[i];
  }

  emxInit_char_T(&c_varargin_2, 2);
  c_varargin_1->data[varargin_1->size[1]] = '\x00';
  i = c_varargin_2->size[0] * c_varargin_2->size[1];
  c_varargin_2->size[0] = 1;
  c_varargin_2->size[1] = varargin_2->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_2, i);
  nbytes = varargin_2->size[1];
  for (i = 0; i < nbytes; i++) {
    c_varargin_2->data[i] = varargin_2->data[i];
  }

  emxInit_char_T(&c_varargin_3, 2);
  c_varargin_2->data[varargin_2->size[1]] = '\x00';
  i = c_varargin_3->size[0] * c_varargin_3->size[1];
  c_varargin_3->size[0] = 1;
  c_varargin_3->size[1] = varargin_3->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_3, i);
  nbytes = varargin_3->size[1];
  for (i = 0; i < nbytes; i++) {
    c_varargin_3->data[i] = varargin_3->data[i];
  }

  emxInit_char_T(&c_varargin_4, 2);
  c_varargin_3->data[varargin_3->size[1]] = '\x00';
  i = c_varargin_4->size[0] * c_varargin_4->size[1];
  c_varargin_4->size[0] = 1;
  c_varargin_4->size[1] = varargin_4->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_4, i);
  nbytes = varargin_4->size[1];
  for (i = 0; i < nbytes; i++) {
    c_varargin_4->data[i] = varargin_4->data[i];
  }

  emxInit_char_T(&c_varargin_5, 2);
  c_varargin_4->data[varargin_4->size[1]] = '\x00';
  i = c_varargin_5->size[0] * c_varargin_5->size[1];
  c_varargin_5->size[0] = 1;
  c_varargin_5->size[1] = varargin_5->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_5, i);
  nbytes = varargin_5->size[1];
  for (i = 0; i < nbytes; i++) {
    c_varargin_5->data[i] = varargin_5->data[i];
  }

  emxInit_char_T(&c_varargin_6, 2);
  c_varargin_5->data[varargin_5->size[1]] = '\x00';
  i = c_varargin_6->size[0] * c_varargin_6->size[1];
  c_varargin_6->size[0] = 1;
  c_varargin_6->size[1] = varargin_6->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_6, i);
  nbytes = varargin_6->size[1];
  for (i = 0; i < nbytes; i++) {
    c_varargin_6->data[i] = varargin_6->data[i];
  }

  emxInit_char_T(&c_varargin_7, 2);
  c_varargin_6->data[varargin_6->size[1]] = '\x00';
  i = c_varargin_7->size[0] * c_varargin_7->size[1];
  c_varargin_7->size[0] = 1;
  c_varargin_7->size[1] = varargin_7->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_7, i);
  nbytes = varargin_7->size[1];
  for (i = 0; i < nbytes; i++) {
    c_varargin_7->data[i] = varargin_7->data[i];
  }

  emxInit_char_T(&c_varargin_8, 2);
  c_varargin_7->data[varargin_7->size[1]] = '\x00';
  i = c_varargin_8->size[0] * c_varargin_8->size[1];
  c_varargin_8->size[0] = 1;
  c_varargin_8->size[1] = varargin_8->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_8, i);
  nbytes = varargin_8->size[1];
  for (i = 0; i < nbytes; i++) {
    c_varargin_8->data[i] = varargin_8->data[i];
  }

  emxInit_char_T(&c_varargin_9, 2);
  c_varargin_8->data[varargin_8->size[1]] = '\x00';
  i = c_varargin_9->size[0] * c_varargin_9->size[1];
  c_varargin_9->size[0] = 1;
  c_varargin_9->size[1] = varargin_9->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_9, i);
  nbytes = varargin_9->size[1];
  for (i = 0; i < nbytes; i++) {
    c_varargin_9->data[i] = varargin_9->data[i];
  }

  emxInit_char_T(&c_varargin_10, 2);
  c_varargin_9->data[varargin_9->size[1]] = '\x00';
  i = c_varargin_10->size[0] * c_varargin_10->size[1];
  c_varargin_10->size[0] = 1;
  c_varargin_10->size[1] = varargin_10->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_10, i);
  nbytes = varargin_10->size[1];
  for (i = 0; i < nbytes; i++) {
    c_varargin_10->data[i] = varargin_10->data[i];
  }

  c_varargin_10->data[varargin_10->size[1]] = '\x00';
  nbytes = snprintf(NULL, 0, "%s%s%s%s%s%s%s%s%s%s", &c_varargin_1->data[0],
                    &c_varargin_2->data[0], &c_varargin_3->data[0],
                    &c_varargin_4->data[0], &c_varargin_5->data[0],
                    &c_varargin_6->data[0], &c_varargin_7->data[0],
                    &c_varargin_8->data[0], &c_varargin_9->data[0],
                    &c_varargin_10->data[0]);
  i = str->size[0] * str->size[1];
  str->size[0] = 1;
  str->size[1] = nbytes + 1;
  emxEnsureCapacity_char_T(str, i);
  snprintf(&str->data[0], (size_t)(nbytes + 1), "%s%s%s%s%s%s%s%s%s%s",
           &b_varargin_1->data[0], &b_varargin_2->data[0], &b_varargin_3->data[0],
           &b_varargin_4->data[0], &b_varargin_5->data[0], &b_varargin_6->data[0],
           &b_varargin_7->data[0], &b_varargin_8->data[0], &b_varargin_9->data[0],
           &b_varargin_10->data[0]);
  i = str->size[0] * str->size[1];
  if (1 > nbytes) {
    str->size[1] = 0;
  } else {
    str->size[1] = nbytes;
  }

  emxEnsureCapacity_char_T(str, i);
  emxFree_char_T(&b_varargin_10);
  emxFree_char_T(&b_varargin_9);
  emxFree_char_T(&b_varargin_8);
  emxFree_char_T(&b_varargin_7);
  emxFree_char_T(&b_varargin_6);
  emxFree_char_T(&b_varargin_5);
  emxFree_char_T(&b_varargin_4);
  emxFree_char_T(&b_varargin_3);
  emxFree_char_T(&b_varargin_2);
  emxFree_char_T(&b_varargin_1);
  emxFree_char_T(&c_varargin_10);
  emxFree_char_T(&c_varargin_9);
  emxFree_char_T(&c_varargin_8);
  emxFree_char_T(&c_varargin_7);
  emxFree_char_T(&c_varargin_6);
  emxFree_char_T(&c_varargin_5);
  emxFree_char_T(&c_varargin_4);
  emxFree_char_T(&c_varargin_3);
  emxFree_char_T(&c_varargin_2);
  emxFree_char_T(&c_varargin_1);
}

void c_sprintf(const emxArray_char_T *varargin_1, const emxArray_char_T
               *varargin_2, emxArray_char_T *str)
{
  emxArray_char_T *b_varargin_1;
  emxArray_char_T *b_varargin_2;
  emxArray_char_T *c_varargin_1;
  emxArray_char_T *c_varargin_2;
  int i;
  int nbytes;
  emxInit_char_T(&b_varargin_1, 2);
  i = b_varargin_1->size[0] * b_varargin_1->size[1];
  b_varargin_1->size[0] = 1;
  b_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_1, i);
  nbytes = varargin_1->size[1];
  for (i = 0; i < nbytes; i++) {
    b_varargin_1->data[i] = varargin_1->data[i];
  }

  emxInit_char_T(&b_varargin_2, 2);
  b_varargin_1->data[varargin_1->size[1]] = '\x00';
  i = b_varargin_2->size[0] * b_varargin_2->size[1];
  b_varargin_2->size[0] = 1;
  b_varargin_2->size[1] = varargin_2->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_2, i);
  nbytes = varargin_2->size[1];
  for (i = 0; i < nbytes; i++) {
    b_varargin_2->data[i] = varargin_2->data[i];
  }

  emxInit_char_T(&c_varargin_1, 2);
  b_varargin_2->data[varargin_2->size[1]] = '\x00';
  i = c_varargin_1->size[0] * c_varargin_1->size[1];
  c_varargin_1->size[0] = 1;
  c_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_1, i);
  nbytes = varargin_1->size[1];
  for (i = 0; i < nbytes; i++) {
    c_varargin_1->data[i] = varargin_1->data[i];
  }

  emxInit_char_T(&c_varargin_2, 2);
  c_varargin_1->data[varargin_1->size[1]] = '\x00';
  i = c_varargin_2->size[0] * c_varargin_2->size[1];
  c_varargin_2->size[0] = 1;
  c_varargin_2->size[1] = varargin_2->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_2, i);
  nbytes = varargin_2->size[1];
  for (i = 0; i < nbytes; i++) {
    c_varargin_2->data[i] = varargin_2->data[i];
  }

  c_varargin_2->data[varargin_2->size[1]] = '\x00';
  nbytes = snprintf(NULL, 0, "%s%s", &c_varargin_1->data[0], &c_varargin_2->
                    data[0]);
  i = str->size[0] * str->size[1];
  str->size[0] = 1;
  str->size[1] = nbytes + 1;
  emxEnsureCapacity_char_T(str, i);
  snprintf(&str->data[0], (size_t)(nbytes + 1), "%s%s", &b_varargin_1->data[0],
           &b_varargin_2->data[0]);
  i = str->size[0] * str->size[1];
  if (1 > nbytes) {
    str->size[1] = 0;
  } else {
    str->size[1] = nbytes;
  }

  emxEnsureCapacity_char_T(str, i);
  emxFree_char_T(&b_varargin_2);
  emxFree_char_T(&b_varargin_1);
  emxFree_char_T(&c_varargin_2);
  emxFree_char_T(&c_varargin_1);
}

/* End of code generation (sprintf.c) */
