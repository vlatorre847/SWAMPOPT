/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * introsort.c
 *
 * Code generation for function 'introsort'
 *
 */

/* Include files */
#include "introsort.h"
#include "get_MILP_types.h"
#include "heapsort.h"
#include "insertionsort.h"
#include "rt_nonfinite.h"

/* Type Definitions */
#ifndef typedef_struct_T
#define typedef_struct_T

typedef struct {
  int xstart;
  int xend;
  int depth;
} struct_T;

#endif                                 /*typedef_struct_T*/

/* Function Definitions */
void introsort(emxArray_int32_T *x, int xend, const cell_wrap_6
               cmp_tunableEnvironment[2])
{
  struct_T st_d_data[120];
  struct_T frame;
  int MAXDEPTH;
  int exitg2;
  int exitg3;
  int i;
  int j;
  int pivot;
  int pmax;
  int pmin;
  int pow2p;
  int st_n;
  boolean_T exitg1;
  boolean_T varargout_1;
  if (1 < xend) {
    if (xend <= 32) {
      insertionsort(x, 1, xend, cmp_tunableEnvironment);
    } else {
      pmax = 31;
      pmin = 0;
      exitg1 = false;
      while ((!exitg1) && (pmax - pmin > 1)) {
        j = (pmin + pmax) >> 1;
        pow2p = 1 << j;
        if (pow2p == xend) {
          pmax = j;
          exitg1 = true;
        } else if (pow2p > xend) {
          pmax = j;
        } else {
          pmin = j;
        }
      }

      MAXDEPTH = (pmax - 1) << 1;
      frame.xstart = 1;
      frame.xend = xend;
      frame.depth = 0;
      pmax = MAXDEPTH << 1;
      for (i = 0; i < pmax; i++) {
        st_d_data[i] = frame;
      }

      st_d_data[0] = frame;
      st_n = 1;
      while (st_n > 0) {
        frame = st_d_data[st_n - 1];
        st_n--;
        i = frame.xend - frame.xstart;
        if (i + 1 <= 32) {
          insertionsort(x, frame.xstart, frame.xend, cmp_tunableEnvironment);
        } else if (frame.depth == MAXDEPTH) {
          b_heapsort(x, frame.xstart, frame.xend, cmp_tunableEnvironment);
        } else {
          pow2p = (frame.xstart + i / 2) - 1;
          i = x->data[frame.xstart - 1];
          pmax = cmp_tunableEnvironment[0].f1->data[x->data[pow2p] - 1];
          pmin = cmp_tunableEnvironment[0].f1->data[i - 1];
          if (pmax < pmin) {
            varargout_1 = true;
          } else if (pmax == pmin) {
            varargout_1 = (cmp_tunableEnvironment[1].f1->data[x->data[pow2p] - 1]
                           < cmp_tunableEnvironment[1].f1->data[i - 1]);
          } else {
            varargout_1 = false;
          }

          if (varargout_1) {
            x->data[frame.xstart - 1] = x->data[pow2p];
            x->data[pow2p] = i;
          }

          if (cmp_tunableEnvironment[0].f1->data[x->data[frame.xend - 1] - 1] <
              cmp_tunableEnvironment[0].f1->data[x->data[frame.xstart - 1] - 1])
          {
            varargout_1 = true;
          } else if (cmp_tunableEnvironment[0].f1->data[x->data[frame.xend - 1]
                     - 1] == cmp_tunableEnvironment[0].f1->data[x->
                     data[frame.xstart - 1] - 1]) {
            varargout_1 = (cmp_tunableEnvironment[1].f1->data[x->data[frame.xend
                           - 1] - 1] < cmp_tunableEnvironment[1].f1->data
                           [x->data[frame.xstart - 1] - 1]);
          } else {
            varargout_1 = false;
          }

          if (varargout_1) {
            pmax = x->data[frame.xstart - 1];
            x->data[frame.xstart - 1] = x->data[frame.xend - 1];
            x->data[frame.xend - 1] = pmax;
          }

          if (cmp_tunableEnvironment[0].f1->data[x->data[frame.xend - 1] - 1] <
              cmp_tunableEnvironment[0].f1->data[x->data[pow2p] - 1]) {
            varargout_1 = true;
          } else if (cmp_tunableEnvironment[0].f1->data[x->data[frame.xend - 1]
                     - 1] == cmp_tunableEnvironment[0].f1->data[x->data[pow2p] -
                     1]) {
            varargout_1 = (cmp_tunableEnvironment[1].f1->data[x->data[frame.xend
                           - 1] - 1] < cmp_tunableEnvironment[1].f1->data
                           [x->data[pow2p] - 1]);
          } else {
            varargout_1 = false;
          }

          if (varargout_1) {
            pmax = x->data[pow2p];
            x->data[pow2p] = x->data[frame.xend - 1];
            x->data[frame.xend - 1] = pmax;
          }

          pivot = x->data[pow2p] - 1;
          x->data[pow2p] = x->data[frame.xend - 2];
          x->data[frame.xend - 2] = pivot + 1;
          pmin = frame.xstart - 1;
          j = frame.xend - 2;
          do {
            exitg2 = 0;
            pmin++;
            do {
              exitg3 = 0;
              i = cmp_tunableEnvironment[0].f1->data[x->data[pmin] - 1];
              if (i < cmp_tunableEnvironment[0].f1->data[pivot]) {
                varargout_1 = true;
              } else if (i == cmp_tunableEnvironment[0].f1->data[pivot]) {
                varargout_1 = (cmp_tunableEnvironment[1].f1->data[x->data[pmin]
                               - 1] < cmp_tunableEnvironment[1].f1->data[pivot]);
              } else {
                varargout_1 = false;
              }

              if (varargout_1) {
                pmin++;
              } else {
                exitg3 = 1;
              }
            } while (exitg3 == 0);

            j--;
            do {
              exitg3 = 0;
              i = cmp_tunableEnvironment[0].f1->data[x->data[j] - 1];
              if (cmp_tunableEnvironment[0].f1->data[pivot] < i) {
                varargout_1 = true;
              } else if (cmp_tunableEnvironment[0].f1->data[pivot] == i) {
                varargout_1 = (cmp_tunableEnvironment[1].f1->data[pivot] <
                               cmp_tunableEnvironment[1].f1->data[x->data[j] - 1]);
              } else {
                varargout_1 = false;
              }

              if (varargout_1) {
                j--;
              } else {
                exitg3 = 1;
              }
            } while (exitg3 == 0);

            if (pmin + 1 >= j + 1) {
              exitg2 = 1;
            } else {
              pmax = x->data[pmin];
              x->data[pmin] = x->data[j];
              x->data[j] = pmax;
            }
          } while (exitg2 == 0);

          x->data[frame.xend - 2] = x->data[pmin];
          x->data[pmin] = pivot + 1;
          if (pmin + 2 < frame.xend) {
            st_d_data[st_n].xstart = pmin + 2;
            st_d_data[st_n].xend = frame.xend;
            st_d_data[st_n].depth = frame.depth + 1;
            st_n++;
          }

          if (frame.xstart < pmin + 1) {
            st_d_data[st_n].xstart = frame.xstart;
            st_d_data[st_n].xend = pmin + 1;
            st_d_data[st_n].depth = frame.depth + 1;
            st_n++;
          }
        }
      }
    }
  }
}

/* End of code generation (introsort.c) */
