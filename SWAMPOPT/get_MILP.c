/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * get_MILP.c
 *
 * Code generation for function 'get_MILP'
 *
 */

/* Include files */
#include "get_MILP.h"
#include "colon.h"
#include "diag.h"
#include "eye.h"
#include "find.h"
#include "get_MILP_data.h"
#include "get_MILP_emxutil.h"
#include "get_MILP_initialize.h"
#include "get_MILP_types.h"
#include "horzcat.h"
#include "kron.h"
#include "minOrMax.h"
#include "nullAssignment.h"
#include "repelem.h"
#include "repmat.h"
#include "repmat1.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "sparse.h"
#include "sparse1.h"
#include "spdiags.h"
#include "sprintf.h"
#include "sum.h"
#include "triu.h"
#include "vertcat.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */
static double rt_powd_snf(double u0, double u1)
{
  double d;
  double d1;
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d == 1.0) {
        y = 1.0;
      } else if (d > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

void get_MILP(int N, int M, int K, int l, int tp_size,
              int tpi_size, int total_delay, const double V0[total_delay], const double R0[l],
              const double H0[l], const int tp[tp_size], const int
              tpi[tpi_size], const double c[l], const double r[N], const
              double b_gamma[l], const double rho[l], const double psi
              [l*l], const double delay[l], double eps_min, const int
              Ki[K], const int Ii[l], const double q[K], const
              int s[K], const int d[K], const double eps[K],
              const double alfa[K], const double beta[K], const
              double Dt[K], const double Dv[K], const double j_funct
              [4], emxArray_real_T *f, emxArray_real_T *A_beg,
              emxArray_real_T *A_nnz, emxArray_int32_T *A_indx,
              emxArray_real_T *A_val, emxArray_char_T *A_sense,
              emxArray_real_T *b, emxArray_real_T *lb, emxArray_real_T *ub,
              emxArray_char_T *ctype, emxArray_real_T *A_j, double *A_r)
{
  d_sparse A102;
  d_sparse A22;
  d_sparse A2eqini;
  d_sparse A3;
  d_sparse A4;
  d_sparse A5;
  d_sparse A6;
  d_sparse A7;
  d_sparse Aeq1;
  d_sparse Aeq5;
  d_sparse Aeq7;
  d_sparse Aeqini;
  d_sparse Aeqini1;
  d_sparse aeq1_temp2;
  d_sparse b_expl_temp;
  d_sparse c_expl_temp;
  d_sparse d_expl_temp;
  d_sparse e_expl_temp;
  d_sparse expl_temp;
  d_sparse f_expl_temp;
  d_sparse g_expl_temp;
  d_sparse h_expl_temp;
  d_sparse i_expl_temp;
  d_sparse j_expl_temp;
  d_sparse k_expl_temp;
  emxArray_boolean_T *c_tempd;
  emxArray_boolean_T *t5_d;
  emxArray_char_T *c_y;
  emxArray_char_T *ctypeE;
  emxArray_char_T *ctypePL;
  emxArray_char_T *ctypeS;
  emxArray_char_T *ctypeZ;
  emxArray_int32_T *appo_tmp;
  emxArray_int32_T *ini_P_tmp;
  emxArray_int32_T *t4_colidx;
  emxArray_int32_T *t4_rowidx;
  emxArray_int32_T *t5_colidx;
  emxArray_int8_T *b21;
  emxArray_int8_T *b22;
  emxArray_int8_T *b4;
  emxArray_int8_T *b8;
  emxArray_int8_T *beq3;
  emxArray_int8_T *beq4;
  emxArray_real_T *Aeq2;
  emxArray_real_T *appo;
  emxArray_real_T *appo1;
  emxArray_real_T *appo_rep;
  emxArray_real_T *b_delay_mat;
  emxArray_real_T *b_r;
  emxArray_real_T *b_temp;
  emxArray_real_T *b_tempd;
  emxArray_real_T *b_templ;
  emxArray_real_T *beq1;
  emxArray_real_T *beqini;
  emxArray_real_T *delay_mat;
  emxArray_real_T *index_canc;
  emxArray_real_T *ini_D;
  emxArray_real_T *ini_L;
  emxArray_real_T *ini_P;
  emxArray_real_T *r1;
  emxArray_real_T *r2;
  emxArray_real_T *t4_d;
  emxArray_real_T *temp;
  emxArray_real_T *temp1;
  emxArray_real_T *tempd;
  emxArray_real_T *tempd1;
  emxArray_real_T *templ;
  emxArray_real_T *templ1;
  double dv[K*K];
  double psi_temp[l*l];
  double b_eps[K];
  double delay_temp[l];
  double delay_data[l-1];
  double a;
  double b_a;
  double b_appo;
  double b_d;
  double b_varargin_1;
  double b_y;
  double c_varargin_1;
  double d1;
  double d2;
  double index_canc_idx_0_tmp;
  double maxval_tmp;
  double varargin_1;
  double varargin_1_tmp;
  double y;
  int b_delay_size[1];
  int delay_size[1];
  int A11_m;
  int A11_n;
  int A12_m;
  int A12_n;
  int Ay_m;
  int Ay_n;
  int aeq1_temp1_m;
  int aeq1_temp1_n;
  int b_i;
  int b_sizes_idx_1;
  int c_i;
  int i;
  int i1;
  int i2;
  int i3;
  int index_canc_idx_0;
  int input_sizes_idx_1;
  int k;
  int loop_ub;
  int loop_ub_tmp;
  int sizes_idx_1;
  int t2_n;
  signed char tmp_data[7];
  boolean_T empty_non_axis_sizes;
  if (!isInitialized_get_MILP) {
    get_MILP_initialize();
  }

  emxInit_real_T(&appo, 2);

  /* N: Number of time intervals */
  /* M: Number of operations */
  /* K: Numer of irrigations */
  /* l: Number of channels */
  /* dt: Time interval duration in minutes */
  /* Ki: Set of the sets of off-takes on the channels */
  /* Ii: Set of the sets of the channels downstream every channel */
  /* q: Quantity of water required by the off-take per time interval */
  /* s: Desidered starting time interval for the irrigation */
  /* d: Desidered duration for the irrigation  */
  /* eps: minimum ratio of water that must be delivered for irrigation */
  /* gamma: seepage ratio for channel */
  /* psi: travelling time from one gate to another */
  /* tp: Time intervals where the gate-keeper cannot operate */
  /* tpi: Time intervals the irrigations cannot start */
  /* c: Maximum inlet volume capacity for every gate */
  /* rho: Minimum ratio of water that must flow through an open gate */
  /* R0: Initial water stored in the channels */
  /* H0: Initial opening of the gates */
  /* tau: delay times for every channel */
  /* delay: time intervals the water takes to go through a channel */
  /* alfa: time priority for the k-th irrigation */
  /* beta: volume priority for the k-th irrigation */
  /* Dt: Maximum delay for an irrigation */
  /* Dv: Minimum volume that can be delivered */
  /* j_funct: Coefficinets of the objective function */
  /* Constraints creation, latex notation used */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  fp1 = fopen("temp.txt", "r"); */
  /*  NN=4*l*N+l^2*N*M+l*N*M+2*K*N; */
  /*  x0_temp=fscanf(fp1,"%f",NN+2*K); */
  /*   V=(vec2mat(x0_temp(1:N*l),l))'; */
  /*   R=(vec2mat(x0_temp(N*l+1:2*N*l),l))'; */
  /*   G=(vec2mat(x0_temp(2*N*l+1:3*N*l),l))'; */
  /*   F=reshape(x0_temp(3*N*l+1:3*N*l+l*l*N*M),l,l,M,N); */
  /*   E=reshape(x0_temp(3*N*l+l*l*N*M+1:3*N*l+l*l*N*M+l*N*M,1),l,M,N); */
  /*   H=reshape(x0_temp(3*N*l+l*l*N*M+l*N*M+K*N+K*N+1:3*N*l+l*l*N*M+l*N*M+K*N+K*N+l*N,1),l,N); */
  /*   */
  /*  fclose(fp1); */
  /* Aeq1 */
  /* (1-\gamma_i)^{\tau_i}V_i^{n-\tau_i} +(1-\gamma_i)R_i^{n-1}- \sum_{k\in K_i} q_kD_k^n-\sum_{j\in I_i}V_j^n-R_i^n=0 */
  i = appo->size[0] * appo->size[1];
  appo->size[0] = (int)l;
  appo->size[1] = (int)l;
  emxEnsureCapacity_real_T(appo, i);
  loop_ub = (int)l * (int)l;
  for (i = 0; i < loop_ub; i++) {
    appo->data[i] = 0.0;
  }

  emxInit_real_T(&appo1, 2);
  i = appo1->size[0] * appo1->size[1];
  appo1->size[0] = (int)l;
  appo1->size[1] = (int)K;
  emxEnsureCapacity_real_T(appo1, i);
  loop_ub = (int)l * (int)K;
  for (i = 0; i < loop_ub; i++) {
    appo1->data[i] = 0.0;
  }

  emxInit_real_T(&temp, 2);
  maxval_tmp = maximum(l,delay);
  i = temp->size[0] * temp->size[1];
  temp->size[0] = (int)l;
  i1 = (int)((maxval_tmp + 1.0) * l);
  temp->size[1] = i1;
  emxEnsureCapacity_real_T(temp, i);
  loop_ub_tmp = (int)l * i1;
  for (i = 0; i < loop_ub_tmp; i++) {
    temp->data[i] = 0.0;
  }

  emxInit_real_T(&templ, 2);
  i = templ->size[0] * templ->size[1];
  templ->size[0] = (int)l;
  templ->size[1] = i1;
  emxEnsureCapacity_real_T(templ, i);
  for (i = 0; i < loop_ub_tmp; i++) {
    templ->data[i] = 0.0;
  }

  emxInit_real_T(&temp1, 2);
  varargin_1_tmp = l * N;
  index_canc_idx_0 = (int)varargin_1_tmp;
  i = temp1->size[0] * temp1->size[1];
  temp1->size[0] = (int)varargin_1_tmp;
  b_d = N + maxval_tmp;
  i1 = (int)(b_d * l);
  temp1->size[1] = i1;
  emxEnsureCapacity_real_T(temp1, i);
  loop_ub_tmp = (int)varargin_1_tmp * i1;
  for (i = 0; i < loop_ub_tmp; i++) {
    temp1->data[i] = 0.0;
  }

  emxInit_real_T(&templ1, 2);
  i = templ1->size[0] * templ1->size[1];
  templ1->size[0] = (int)varargin_1_tmp;
  templ1->size[1] = i1;
  emxEnsureCapacity_real_T(templ1, i);
  for (i = 0; i < loop_ub_tmp; i++) {
    templ1->data[i] = 0.0;
  }

  emxInit_real_T(&tempd, 2);
  i = tempd->size[0] * tempd->size[1];
  tempd->size[0] = (int)l;
  i1 = (int)(K * (maxval_tmp + 1.0));
  tempd->size[1] = i1;
  emxEnsureCapacity_real_T(tempd, i);
  loop_ub = (int)l * i1;
  for (i = 0; i < loop_ub; i++) {
    tempd->data[i] = 0.0;
  }

  emxInit_real_T(&tempd1, 2);
  i = tempd1->size[0] * tempd1->size[1];
  tempd1->size[0] = (int)varargin_1_tmp;
  i1 = (int)(K * b_d);
  tempd1->size[1] = i1;
  emxEnsureCapacity_real_T(tempd1, i);
  loop_ub = (int)varargin_1_tmp * i1;
  for (i = 0; i < loop_ub; i++) {
    tempd1->data[i] = 0.0;
  }

  i = (int)l;
  for (b_i = 0; b_i < i; b_i++) {
    k = 0;
    sizes_idx_1 = 0;
    for (c_i = 0; c_i < l; c_i++) {
      if (Ii[c_i] == (double)b_i + 1.0) {
        k++;
        tmp_data[sizes_idx_1] = (signed char)(c_i + 1);
        sizes_idx_1++;
      }
    }

    for (i1 = 0; i1 < k; i1++) {
      appo->data[b_i + appo->size[0] * (tmp_data[i1] - 1)] = -1.0;
    }

    i1 = (int)K;
    for (sizes_idx_1 = 0; sizes_idx_1 < i1; sizes_idx_1++) {
      if (Ki[sizes_idx_1] == (double)b_i + 1.0) {
        appo1->data[b_i + appo1->size[0] * sizes_idx_1] = -q[sizes_idx_1];
      }
    }

    b_d = delay[b_i];
    d1 = b_d * l;
    if (d1 + 1.0 > (b_d + 1.0) * l) {
      i1 = 1;
    } else {
      i1 = (int)(d1 + 1.0);
    }

    loop_ub = appo->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      temp->data[b_i + temp->size[0] * ((i1 + i2) - 1)] = appo->data[b_i +
        appo->size[0] * i2];
    }

    templ->data[b_i + templ->size[0] * ((int)(d1 + ((double)b_i + 1.0)) - 1)] =
      1.0;
    b_d = delay[b_i];
    d1 = b_d * K + 1.0;
    if (d1 > (b_d + 1.0) * K) {
      i1 = 1;
    } else {
      i1 = (int)d1;
    }

    loop_ub = appo1->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      tempd->data[b_i + tempd->size[0] * ((i1 + i2) - 1)] = appo1->data[b_i +
        appo1->size[0] * i2];
    }
  }

  i = (int)N;
  if (0 <= (int)N - 1) {
    for (i1 = 0; i1 < l; i1++) {
      delay_temp[i1] = delay[i1] + 1.0;
    }
  }

  emxInit_real_T(&appo_rep, 1);
  emxInit_int32_T(&ini_P_tmp, 1);
  for (input_sizes_idx_1 = 0; input_sizes_idx_1 < i; input_sizes_idx_1++) {
    b_appo = (((double)input_sizes_idx_1 + 1.0) - 1.0) * l;
    b_d = ((double)input_sizes_idx_1 + 1.0) * l;
    if (b_appo + 1.0 > b_d) {
      i1 = 1;
    } else {
      i1 = (int)(b_appo + 1.0);
    }

    d1 = maximum(l,delay_temp);
    loop_ub = (int)floor(d1 * l - 1.0);
    i2 = appo_rep->size[0];
    appo_rep->size[0] = loop_ub + 1;
    emxEnsureCapacity_real_T(appo_rep, i2);
    for (i2 = 0; i2 <= loop_ub; i2++) {
      appo_rep->data[i2] = b_appo + ((double)i2 + 1.0);
    }

    loop_ub = temp->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      sizes_idx_1 = temp->size[0];
      for (Ay_n = 0; Ay_n < sizes_idx_1; Ay_n++) {
        temp1->data[((i1 + Ay_n) + temp1->size[0] * ((int)appo_rep->data[i2] - 1))
          - 1] = temp->data[Ay_n + temp->size[0] * i2];
      }
    }

    y = (((double)input_sizes_idx_1 + 1.0) - 1.0) * l + 1.0;
    if (y > b_d) {
      i1 = 1;
    } else {
      i1 = (int)y;
    }

    loop_ub = templ->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      sizes_idx_1 = templ->size[0];
      for (Ay_n = 0; Ay_n < sizes_idx_1; Ay_n++) {
        templ1->data[((i1 + Ay_n) + templ1->size[0] * ((int)appo_rep->data[i2] -
          1)) - 1] = templ->data[Ay_n + templ->size[0] * i2];
      }
    }

    y = (((double)input_sizes_idx_1 + 1.0) - 1.0) * l + 1.0;
    if (y > b_d) {
      i1 = 1;
    } else {
      i1 = (int)y;
    }

    b_appo = (((double)input_sizes_idx_1 + 1.0) - 1.0) * K;
    loop_ub = (int)floor(d1 * K - 1.0);
    i2 = ini_P_tmp->size[0];
    ini_P_tmp->size[0] = loop_ub + 1;
    emxEnsureCapacity_int32_T(ini_P_tmp, i2);
    for (i2 = 0; i2 <= loop_ub; i2++) {
      ini_P_tmp->data[i2] = (int)(b_appo + ((double)i2 + 1.0)) - 1;
    }

    loop_ub = tempd->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      sizes_idx_1 = tempd->size[0];
      for (Ay_n = 0; Ay_n < sizes_idx_1; Ay_n++) {
        tempd1->data[((i1 + Ay_n) + tempd1->size[0] * ini_P_tmp->data[i2]) - 1] =
          tempd->data[Ay_n + tempd->size[0] * i2];
      }
    }
  }

  emxInit_int32_T(&appo_tmp, 2);
  i = appo_tmp->size[0] * appo_tmp->size[1];
  appo_tmp->size[0] = 1;
  appo_tmp->size[1] = (int)floor((double)temp1->size[1] - (varargin_1_tmp + 1.0))
    + 1;
  emxEnsureCapacity_int32_T(appo_tmp, i);
  loop_ub = (int)floor((double)temp1->size[1] - (varargin_1_tmp + 1.0));
  for (i = 0; i <= loop_ub; i++) {
    appo_tmp->data[i] = (int)((varargin_1_tmp + 1.0) + (double)i);
  }

  b_nullAssignment(temp1, appo_tmp);
  b_d = N * l + 1.0;
  i = appo_tmp->size[0] * appo_tmp->size[1];
  appo_tmp->size[0] = 1;
  appo_tmp->size[1] = (int)floor((double)templ1->size[1] - b_d) + 1;
  emxEnsureCapacity_int32_T(appo_tmp, i);
  loop_ub = (int)floor((double)templ1->size[1] - b_d);
  for (i = 0; i <= loop_ub; i++) {
    appo_tmp->data[i] = (int)(b_d + (double)i);
  }

  b_nullAssignment(templ1, appo_tmp);
  varargin_1 = K * N;
  i = appo_tmp->size[0] * appo_tmp->size[1];
  appo_tmp->size[0] = 1;
  appo_tmp->size[1] = (int)floor((double)tempd1->size[1] - (varargin_1 + 1.0)) +
    1;
  emxEnsureCapacity_int32_T(appo_tmp, i);
  loop_ub = (int)floor((double)tempd1->size[1] - (varargin_1 + 1.0));
  for (i = 0; i <= loop_ub; i++) {
    appo_tmp->data[i] = (int)((varargin_1 + 1.0) + (double)i);
  }

  emxInit_real_T(&ini_P, 2);
  b_nullAssignment(tempd1, appo_tmp);
  i = (int)(l * (maxval_tmp - 1.0));
  i1 = ini_P->size[0] * ini_P->size[1];
  ini_P->size[0] = i;
  ini_P->size[1] = i;
  emxEnsureCapacity_real_T(ini_P, i1);
  loop_ub = i * i;
  for (i1 = 0; i1 < loop_ub; i1++) {
    ini_P->data[i1] = 0.0;
  }

  emxInit_real_T(&ini_D, 2);
  i1 = ini_D->size[0] * ini_D->size[1];
  ini_D->size[0] = i;
  i2 = (int)(K * (maxval_tmp - 1.0));
  ini_D->size[1] = i2;
  emxEnsureCapacity_real_T(ini_D, i1);
  loop_ub = i * i2;
  for (i1 = 0; i1 < loop_ub; i1++) {
    ini_D->data[i1] = 0.0;
  }

  emxInit_real_T(&ini_L, 2);
  i1 = ini_L->size[0] * ini_L->size[1];
  ini_L->size[0] = i;
  Ay_n = (int)(l * maxval_tmp);
  ini_L->size[1] = Ay_n;
  emxEnsureCapacity_real_T(ini_L, i1);
  loop_ub = i * Ay_n;
  for (i = 0; i < loop_ub; i++) {
    ini_L->data[i] = 0.0;
  }

  eye(l, templ);
  i = tempd->size[0] * tempd->size[1];
  tempd->size[0] = (int)l;
  tempd->size[1] = (int)(maxval_tmp - 1.0);
  emxEnsureCapacity_real_T(tempd, i);
  loop_ub = (int)l * (int)(maxval_tmp - 1.0);
  for (i = 0; i < loop_ub; i++) {
    tempd->data[i] = 0.0;
  }

  emxInit_real_T(&delay_mat, 3);
  i = delay_mat->size[0] * delay_mat->size[1] * delay_mat->size[2];
  delay_mat->size[0] = (int)(maxval_tmp - 1.0);
  delay_mat->size[1] = (int)(maxval_tmp - 1.0);
  delay_mat->size[2] = (int)l;
  emxEnsureCapacity_real_T(delay_mat, i);
  loop_ub = (int)(maxval_tmp - 1.0) * (int)(maxval_tmp - 1.0) * (int)l;
  for (i = 0; i < loop_ub; i++) {
    delay_mat->data[i] = 0.0;
  }

  i = (int)l;
  emxInit_real_T(&Aeq2, 2);
  emxInit_real_T(&b_tempd, 2);
  emxInit_real_T(&b_r, 2);
  emxInit_real_T(&b_delay_mat, 2);
  for (b_i = 0; b_i < i; b_i++) {
    b_d = delay[b_i];
    if (1.0 > b_d - 1.0) {
      loop_ub = 0;
    } else {
      loop_ub = (int)(b_d - 1.0);
    }

    for (i1 = 0; i1 < loop_ub; i1++) {
      tempd->data[b_i + tempd->size[0] * i1] = 1.0;
    }

    loop_ub = tempd->size[1];
    i1 = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = tempd->size[1];
    emxEnsureCapacity_real_T(b_tempd, i1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_tempd->data[i1] = tempd->data[b_i + tempd->size[0] * i1];
    }

    diag(b_tempd, b_r);
    loop_ub = b_r->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      sizes_idx_1 = b_r->size[0];
      for (Ay_n = 0; Ay_n < sizes_idx_1; Ay_n++) {
        delay_mat->data[(Ay_n + delay_mat->size[0] * i1) + delay_mat->size[0] *
          delay_mat->size[1] * b_i] = b_r->data[Ay_n + b_r->size[0] * i1];
      }
    }

    b_d = (((double)b_i + 1.0) - 1.0) * (maxval_tmp - 1.0) + 1.0;
    d1 = ((double)b_i + 1.0) * (maxval_tmp - 1.0);
    if (b_d > d1) {
      i1 = 1;
    } else {
      i1 = (int)b_d;
    }

    loop_ub = delay_mat->size[0];
    sizes_idx_1 = delay_mat->size[1];
    Ay_n = b_delay_mat->size[0] * b_delay_mat->size[1];
    b_delay_mat->size[0] = delay_mat->size[0];
    b_delay_mat->size[1] = delay_mat->size[1];
    emxEnsureCapacity_real_T(b_delay_mat, Ay_n);
    for (Ay_n = 0; Ay_n < sizes_idx_1; Ay_n++) {
      for (Ay_m = 0; Ay_m < loop_ub; Ay_m++) {
        b_delay_mat->data[Ay_m + b_delay_mat->size[0] * Ay_n] = delay_mat->data
          [(Ay_m + delay_mat->size[0] * Ay_n) + delay_mat->size[0] *
          delay_mat->size[1] * b_i];
      }
    }

    loop_ub = appo->size[1];
    Ay_n = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = appo->size[1];
    emxEnsureCapacity_real_T(b_tempd, Ay_n);
    for (Ay_n = 0; Ay_n < loop_ub; Ay_n++) {
      b_tempd->data[Ay_n] = appo->data[b_i + appo->size[0] * Ay_n];
    }

    kron(b_delay_mat, b_tempd, b_r);
    loop_ub = b_r->size[1];
    for (Ay_n = 0; Ay_n < loop_ub; Ay_n++) {
      sizes_idx_1 = b_r->size[0];
      for (Ay_m = 0; Ay_m < sizes_idx_1; Ay_m++) {
        ini_P->data[((i1 + Ay_m) + ini_P->size[0] * Ay_n) - 1] = b_r->data[Ay_m
          + b_r->size[0] * Ay_n];
      }
    }

    if (b_d > d1) {
      i1 = 1;
    } else {
      i1 = (int)b_d;
    }

    loop_ub = delay_mat->size[0];
    sizes_idx_1 = delay_mat->size[1];
    Ay_n = b_delay_mat->size[0] * b_delay_mat->size[1];
    b_delay_mat->size[0] = delay_mat->size[0];
    b_delay_mat->size[1] = delay_mat->size[1];
    emxEnsureCapacity_real_T(b_delay_mat, Ay_n);
    for (Ay_n = 0; Ay_n < sizes_idx_1; Ay_n++) {
      for (Ay_m = 0; Ay_m < loop_ub; Ay_m++) {
        b_delay_mat->data[Ay_m + b_delay_mat->size[0] * Ay_n] = delay_mat->data
          [(Ay_m + delay_mat->size[0] * Ay_n) + delay_mat->size[0] *
          delay_mat->size[1] * b_i];
      }
    }

    loop_ub = appo1->size[1];
    Ay_n = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = appo1->size[1];
    emxEnsureCapacity_real_T(b_tempd, Ay_n);
    for (Ay_n = 0; Ay_n < loop_ub; Ay_n++) {
      b_tempd->data[Ay_n] = appo1->data[b_i + appo1->size[0] * Ay_n];
    }

    kron(b_delay_mat, b_tempd, b_r);
    loop_ub = b_r->size[1];
    for (Ay_n = 0; Ay_n < loop_ub; Ay_n++) {
      sizes_idx_1 = b_r->size[0];
      for (Ay_m = 0; Ay_m < sizes_idx_1; Ay_m++) {
        ini_D->data[((i1 + Ay_m) + ini_D->size[0] * Ay_n) - 1] = b_r->data[Ay_m
          + b_r->size[0] * Ay_n];
      }
    }

    loop_ub = delay_mat->size[0];
    sizes_idx_1 = delay_mat->size[1];
    i1 = b_delay_mat->size[0] * b_delay_mat->size[1];
    b_delay_mat->size[0] = delay_mat->size[0];
    b_delay_mat->size[1] = delay_mat->size[1];
    emxEnsureCapacity_real_T(b_delay_mat, i1);
    for (i1 = 0; i1 < sizes_idx_1; i1++) {
      for (Ay_n = 0; Ay_n < loop_ub; Ay_n++) {
        b_delay_mat->data[Ay_n + b_delay_mat->size[0] * i1] = delay_mat->data
          [(Ay_n + delay_mat->size[0] * i1) + delay_mat->size[0] *
          delay_mat->size[1] * b_i] * (1.0 - b_gamma[b_i]);
      }
    }

    loop_ub = templ->size[1];
    i1 = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = templ->size[1];
    emxEnsureCapacity_real_T(b_tempd, i1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_tempd->data[i1] = templ->data[b_i + templ->size[0] * i1];
    }

    kron(b_delay_mat, b_tempd, Aeq2);
    if ((Aeq2->size[0] != 0) && (Aeq2->size[1] != 0)) {
      c_i = Aeq2->size[0];
    } else if (((int)(maxval_tmp - 1.0) != 0) && ((int)l != 0)) {
      c_i = (int)(maxval_tmp - 1.0);
    } else {
      c_i = Aeq2->size[0];
      if (c_i <= 0) {
        c_i = 0;
      }

      if ((int)(maxval_tmp - 1.0) > c_i) {
        c_i = (int)(maxval_tmp - 1.0);
      }
    }

    empty_non_axis_sizes = (c_i == 0);
    if (empty_non_axis_sizes || ((Aeq2->size[0] != 0) && (Aeq2->size[1] != 0)))
    {
      input_sizes_idx_1 = Aeq2->size[1];
    } else {
      input_sizes_idx_1 = 0;
    }

    if (empty_non_axis_sizes || (((int)(maxval_tmp - 1.0) != 0) && ((int)l != 0)))
    {
      b_sizes_idx_1 = (int)l;
    } else {
      b_sizes_idx_1 = 0;
    }

    loop_ub = delay_mat->size[0];
    sizes_idx_1 = delay_mat->size[1];
    i1 = b_delay_mat->size[0] * b_delay_mat->size[1];
    b_delay_mat->size[0] = delay_mat->size[0];
    b_delay_mat->size[1] = delay_mat->size[1];
    emxEnsureCapacity_real_T(b_delay_mat, i1);
    for (i1 = 0; i1 < sizes_idx_1; i1++) {
      for (Ay_n = 0; Ay_n < loop_ub; Ay_n++) {
        b_delay_mat->data[Ay_n + b_delay_mat->size[0] * i1] = delay_mat->data
          [(Ay_n + delay_mat->size[0] * i1) + delay_mat->size[0] *
          delay_mat->size[1] * b_i];
      }
    }

    loop_ub = templ->size[1];
    i1 = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = templ->size[1];
    emxEnsureCapacity_real_T(b_tempd, i1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_tempd->data[i1] = templ->data[b_i + templ->size[0] * i1];
    }

    kron(b_delay_mat, b_tempd, temp);
    if (((int)(maxval_tmp - 1.0) != 0) && ((int)l != 0)) {
      k = (int)(maxval_tmp - 1.0);
    } else if ((temp->size[0] != 0) && (temp->size[1] != 0)) {
      k = temp->size[0];
    } else {
      k = (int)(maxval_tmp - 1.0);
      if (k <= 0) {
        k = 0;
      }

      if (temp->size[0] > k) {
        k = temp->size[0];
      }
    }

    empty_non_axis_sizes = (k == 0);
    if (empty_non_axis_sizes || (((int)(maxval_tmp - 1.0) != 0) && ((int)l != 0)))
    {
      loop_ub = (int)l;
    } else {
      loop_ub = 0;
    }

    if (empty_non_axis_sizes || ((temp->size[0] != 0) && (temp->size[1] != 0)))
    {
      sizes_idx_1 = temp->size[1];
    } else {
      sizes_idx_1 = 0;
    }

    if (b_d > d1) {
      i1 = 1;
    } else {
      i1 = (int)b_d;
    }

    Ay_n = b_delay_mat->size[0] * b_delay_mat->size[1];
    b_delay_mat->size[0] = c_i;
    b_delay_mat->size[1] = input_sizes_idx_1 + b_sizes_idx_1;
    emxEnsureCapacity_real_T(b_delay_mat, Ay_n);
    for (Ay_n = 0; Ay_n < input_sizes_idx_1; Ay_n++) {
      for (Ay_m = 0; Ay_m < c_i; Ay_m++) {
        b_delay_mat->data[Ay_m + b_delay_mat->size[0] * Ay_n] = Aeq2->data[Ay_m
          + c_i * Ay_n];
      }
    }

    for (Ay_n = 0; Ay_n < b_sizes_idx_1; Ay_n++) {
      for (Ay_m = 0; Ay_m < c_i; Ay_m++) {
        b_delay_mat->data[Ay_m + b_delay_mat->size[0] * (Ay_n +
          input_sizes_idx_1)] = 0.0;
      }
    }

    Ay_n = b_r->size[0] * b_r->size[1];
    b_r->size[0] = k;
    b_r->size[1] = loop_ub + sizes_idx_1;
    emxEnsureCapacity_real_T(b_r, Ay_n);
    for (Ay_n = 0; Ay_n < loop_ub; Ay_n++) {
      for (Ay_m = 0; Ay_m < k; Ay_m++) {
        b_r->data[Ay_m + b_r->size[0] * Ay_n] = 0.0;
      }
    }

    for (Ay_n = 0; Ay_n < sizes_idx_1; Ay_n++) {
      for (Ay_m = 0; Ay_m < k; Ay_m++) {
        b_r->data[Ay_m + b_r->size[0] * (Ay_n + loop_ub)] = temp->data[Ay_m + k *
          Ay_n];
      }
    }

    loop_ub = b_delay_mat->size[1];
    for (Ay_n = 0; Ay_n < loop_ub; Ay_n++) {
      sizes_idx_1 = b_delay_mat->size[0];
      for (Ay_m = 0; Ay_m < sizes_idx_1; Ay_m++) {
        ini_L->data[((i1 + Ay_m) + ini_L->size[0] * Ay_n) - 1] =
          b_delay_mat->data[Ay_m + b_delay_mat->size[0] * Ay_n] - b_r->data[Ay_m
          + b_r->size[0] * Ay_n];
      }
    }
  }

  emxFree_real_T(&delay_mat);
  emxInit_boolean_T(&c_tempd, 2);
  i = c_tempd->size[0] * c_tempd->size[1];
  c_tempd->size[0] = tempd->size[1];
  c_tempd->size[1] = tempd->size[0];
  emxEnsureCapacity_boolean_T(c_tempd, i);
  loop_ub = tempd->size[0];
  for (i = 0; i < loop_ub; i++) {
    sizes_idx_1 = tempd->size[1];
    for (i1 = 0; i1 < sizes_idx_1; i1++) {
      c_tempd->data[i1 + c_tempd->size[0] * i] = (tempd->data[i + tempd->size[0]
        * i1] == 0.0);
    }
  }

  eml_find(c_tempd, ini_P_tmp);
  nullAssignment(ini_P, ini_P_tmp);
  nullAssignment(ini_D, ini_P_tmp);
  nullAssignment(ini_L, ini_P_tmp);
  emxFree_boolean_T(&c_tempd);
  for (k = 0; k < l; k++) {
    delay_temp[k] = rt_powd_snf(1.0 - b_gamma[k], delay[k]);
  }

  b_diag(l,delay_temp, psi_temp);
  emxInit_real_T(&r1, 2);
  eye(N, r1);
  b_kron(l,r1, psi_temp, b_r);
  loop_ub = temp1->size[0] * temp1->size[1];
  for (i = 0; i < loop_ub; i++) {
    temp1->data[i] += b_r->data[i];
  }

  emxInitStruct_sparse(&A102);
  emxInitStruct_sparse(&expl_temp);
  emxInitStruct_sparse(&b_expl_temp);
  emxInitStruct_sparse(&c_expl_temp);
  emxInitStruct_sparse(&d_expl_temp);
  emxInitStruct_sparse(&e_expl_temp);
  emxInit_real_T(&b_temp, 2);

  /* Aeqini and Aeqini1 */
  /*  \sum_{k\in K_i} q_kD_k^n+\sum_{j\in I_i}V_j^n+R_i^n= (1-\gamma_i)R^n_i \forall \, i=1,\dots,l, \forall \, n=0,\dots,\max\{\tau_i\}, */
  sparse(l, l * (N - 1.0), e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx,
         &b_i, &t2_n, &c_i);
  eye(l, temp);
  sparse(l, l * (N - 1.0), d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  sparse(l, (N * l + l * l * N * M) + l * N * M, c_expl_temp.d,
         c_expl_temp.colidx, c_expl_temp.rowidx, &Ay_m, &Ay_n, &c_i);
  sparse(l, N * K, b_expl_temp.d, b_expl_temp.colidx, b_expl_temp.rowidx, &k,
         &sizes_idx_1, &c_i);
  sparse(l, (N - 1.0) * K, A102.d, A102.colidx, A102.rowidx, &A102.m, &A102.n,
         &A102.maxnz);
  sparse(l, N * l, expl_temp.d, expl_temp.colidx, expl_temp.rowidx,
         &aeq1_temp1_m, &aeq1_temp1_n, &c_i);
  i = b_temp->size[0] * b_temp->size[1];
  b_temp->size[0] = temp->size[0];
  b_temp->size[1] = temp->size[1];
  emxEnsureCapacity_real_T(b_temp, i);
  loop_ub = temp->size[0] * temp->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_temp->data[i] = -temp->data[i];
  }

  emxInitStruct_sparse(&Aeqini);
  emxInitStruct_sparse(&Aeqini1);
  emxInitStruct_sparse(&aeq1_temp2);
  sparse_horzcat(appo, b_i, t2_n, b_temp, input_sizes_idx_1, b_sizes_idx_1, Ay_m,
                 Ay_n, k, sizes_idx_1, appo1, A102.m, A102.n, aeq1_temp1_m,
                 aeq1_temp1_n, &Aeqini);
  sparse(ini_P->size[0], l, e_expl_temp.d, e_expl_temp.colidx,
         e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  sparse(ini_P->size[0], l * (N - maxval_tmp), d_expl_temp.d, d_expl_temp.colidx,
         d_expl_temp.rowidx, &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  sparse(ini_P->size[0], (N - maxval_tmp) * l, c_expl_temp.d, c_expl_temp.colidx,
         c_expl_temp.rowidx, &Ay_m, &Ay_n, &c_i);
  sparse(ini_P->size[0], ((N * l + l * l * N * M) + l * N * M) + K * N,
         b_expl_temp.d, b_expl_temp.colidx, b_expl_temp.rowidx, &k, &sizes_idx_1,
         &c_i);
  sparse(ini_P->size[0], K, A102.d, A102.colidx, A102.rowidx, &A102.m, &A102.n,
         &A102.maxnz);
  sparse(ini_P->size[0], K * (N - maxval_tmp), expl_temp.d, expl_temp.colidx,
         expl_temp.rowidx, &aeq1_temp1_m, &aeq1_temp1_n, &c_i);
  sparse(ini_P->size[0], N * l, aeq1_temp2.d, aeq1_temp2.colidx,
         aeq1_temp2.rowidx, &aeq1_temp2.m, &aeq1_temp2.n, &aeq1_temp2.maxnz);
  b_sparse_horzcat(b_i, t2_n, ini_P, input_sizes_idx_1, b_sizes_idx_1, ini_L,
                   Ay_m, Ay_n, k, sizes_idx_1, A102.m, A102.n, ini_D,
                   aeq1_temp1_m, aeq1_temp1_n, aeq1_temp2, &Aeqini1);
  emxFree_real_T(&ini_L);
  emxFree_real_T(&ini_D);
  emxFree_real_T(&appo1);
  if (l + 1.0 > varargin_1_tmp) {
    i = -1;
    i1 = -1;
  } else {
    i = (int)(l + 1.0) - 2;
    i1 = (int)varargin_1_tmp - 1;
  }

  for (Ay_n = 0; Ay_n < l; Ay_n++) {
    delay_temp[Ay_n] = 1.0 - b_gamma[Ay_n];
  }

  repmat(l,delay_temp, N, appo_rep);
  loop_ub = i1 - i;
  if ((templ1->size[0] != 0) && (loop_ub != 0)) {
    k = templ1->size[0];
  } else if (((int)varargin_1_tmp != 0) && ((int)l != 0)) {
    k = (int)varargin_1_tmp;
  } else {
    if (templ1->size[0] > 0) {
      k = templ1->size[0];
    } else {
      k = 0;
    }

    if ((int)varargin_1_tmp > k) {
      k = (int)varargin_1_tmp;
    }
  }

  empty_non_axis_sizes = (k == 0);
  if (empty_non_axis_sizes || ((templ1->size[0] != 0) && (loop_ub != 0))) {
    input_sizes_idx_1 = i1 - i;
  } else {
    input_sizes_idx_1 = 0;
  }

  if (empty_non_axis_sizes || (((int)varargin_1_tmp != 0) && ((int)l != 0))) {
    b_sizes_idx_1 = (int)l;
  } else {
    b_sizes_idx_1 = 0;
  }

  sizes_idx_1 = templ1->size[0];
  i1 = b_delay_mat->size[0] * b_delay_mat->size[1];
  b_delay_mat->size[0] = templ1->size[0];
  b_delay_mat->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_delay_mat, i1);
  for (i1 = 0; i1 < loop_ub; i1++) {
    for (Ay_n = 0; Ay_n < sizes_idx_1; Ay_n++) {
      b_delay_mat->data[Ay_n + b_delay_mat->size[0] * i1] = templ1->data[Ay_n +
        templ1->size[0] * ((i + i1) + 1)];
    }
  }

  i = templ->size[0] * templ->size[1];
  templ->size[0] = k;
  templ->size[1] = input_sizes_idx_1 + b_sizes_idx_1;
  emxEnsureCapacity_real_T(templ, i);
  for (i = 0; i < input_sizes_idx_1; i++) {
    for (i1 = 0; i1 < k; i1++) {
      templ->data[i1 + templ->size[0] * i] = b_delay_mat->data[i1 + k * i];
    }
  }

  for (i = 0; i < b_sizes_idx_1; i++) {
    for (i1 = 0; i1 < k; i1++) {
      templ->data[i1 + templ->size[0] * (i + input_sizes_idx_1)] = 0.0;
    }
  }

  emxInit_real_T(&b_templ, 1);
  for (b_i = 0; b_i < index_canc_idx_0; b_i++) {
    sizes_idx_1 = templ->size[0] - 1;
    i = b_templ->size[0];
    b_templ->size[0] = templ->size[0];
    emxEnsureCapacity_real_T(b_templ, i);
    for (i = 0; i <= sizes_idx_1; i++) {
      b_templ->data[i] = templ->data[i + templ->size[0] * b_i] * appo_rep->
        data[i];
    }

    loop_ub = b_templ->size[0];
    for (i = 0; i < loop_ub; i++) {
      templ->data[i + templ->size[0] * b_i] = b_templ->data[i];
    }
  }

  sparse(N * l, N * l, e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx,
         &b_i, &t2_n, &c_i);
  sparse(N * l, l * l * N * M, d_expl_temp.d, d_expl_temp.colidx,
         d_expl_temp.rowidx, &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  sparse(N * l, l * N * M, c_expl_temp.d, c_expl_temp.colidx, c_expl_temp.rowidx,
         &Ay_m, &Ay_n, &c_i);
  sparse(N * l, K * N, b_expl_temp.d, b_expl_temp.colidx, b_expl_temp.rowidx, &k,
         &sizes_idx_1, &c_i);
  sparse(N * l, N * l, A102.d, A102.colidx, A102.rowidx, &A102.m, &A102.n,
         &A102.maxnz);
  i = b_delay_mat->size[0] * b_delay_mat->size[1];
  b_delay_mat->size[0] = templ->size[0];
  b_delay_mat->size[1] = templ->size[1];
  emxEnsureCapacity_real_T(b_delay_mat, i);
  loop_ub = templ->size[0] * templ->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_delay_mat->data[i] = templ->data[i] - templ1->data[i];
  }

  emxFree_real_T(&templ1);
  emxInitStruct_sparse(&Aeq1);
  c_sparse_horzcat(temp1, b_delay_mat, b_i, t2_n, input_sizes_idx_1,
                   b_sizes_idx_1, Ay_m, Ay_n, k, sizes_idx_1, tempd1, A102.m,
                   A102.n, &Aeq1);
  emxFree_real_T(&tempd1);
  emxFree_real_T(&temp1);
  for (b_i = 0; b_i < l; b_i++) {
    delay_temp[b_i] = N - delay[b_i];
  }

  emxInit_real_T(&index_canc, 2);
  index_canc->size[0] = 0;
  index_canc->size[1] = 1;
  i = (int)l;
  for (b_i = 0; b_i < i; b_i++) {
    b_d = delay_temp[b_i];
    if (rtIsNaN(b_d) || rtIsNaN(N - 1.0)) {
      i1 = b_tempd->size[0] * b_tempd->size[1];
      b_tempd->size[0] = 1;
      b_tempd->size[1] = 1;
      emxEnsureCapacity_real_T(b_tempd, i1);
      b_tempd->data[0] = rtNaN;
    } else if (N - 1.0 < b_d) {
      b_tempd->size[0] = 1;
      b_tempd->size[1] = 0;
    } else if ((rtIsInf(b_d) || rtIsInf(N - 1.0)) && (b_d == N - 1.0)) {
      i1 = b_tempd->size[0] * b_tempd->size[1];
      b_tempd->size[0] = 1;
      b_tempd->size[1] = 1;
      emxEnsureCapacity_real_T(b_tempd, i1);
      b_tempd->data[0] = rtNaN;
    } else if (floor(b_d) == b_d) {
      i1 = b_tempd->size[0] * b_tempd->size[1];
      b_tempd->size[0] = 1;
      loop_ub = (int)floor((N - 1.0) - b_d);
      b_tempd->size[1] = loop_ub + 1;
      emxEnsureCapacity_real_T(b_tempd, i1);
      for (i1 = 0; i1 <= loop_ub; i1++) {
        b_tempd->data[i1] = delay_temp[b_i] + (double)i1;
      }
    } else {
      eml_float_colon(b_d, N - 1.0, b_tempd);
    }

    i1 = appo_rep->size[0];
    appo_rep->size[0] = index_canc->size[0] + b_tempd->size[1];
    emxEnsureCapacity_real_T(appo_rep, i1);
    loop_ub = index_canc->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      appo_rep->data[i1] = index_canc->data[i1];
    }

    loop_ub = b_tempd->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      appo_rep->data[i1 + index_canc->size[0]] = b_tempd->data[i1] * l +
        ((double)b_i + 1.0);
    }

    i1 = index_canc->size[0] * index_canc->size[1];
    index_canc->size[0] = appo_rep->size[0];
    index_canc->size[1] = 1;
    emxEnsureCapacity_real_T(index_canc, i1);
    loop_ub = appo_rep->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      index_canc->data[i1] = appo_rep->data[i1];
    }
  }

  i = appo_rep->size[0];
  appo_rep->size[0] = index_canc->size[0];
  emxEnsureCapacity_real_T(appo_rep, i);
  loop_ub = index_canc->size[0];
  for (i = 0; i < loop_ub; i++) {
    appo_rep->data[i] = index_canc->data[i];
  }

  emxFree_real_T(&index_canc);
  sort(appo_rep);

  /* VVVVVVVVVVVVVVVVVV  Changes from here   VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV */
  b_d = delay[0];
  for (k = 0; k < 6; k++) {
    b_d += delay[k + 1];
  }

  i = (int)(((-1.0 - b_d) + 1.0) / -1.0);
  for (b_i = 0; b_i < i; b_i++) {
    b_appo = b_d + -(double)b_i;
    d1 = appo_rep->data[(int)b_appo - 1];
    if (rtIsNaN(d1 - 1.0)) {
      i1 = b_tempd->size[0] * b_tempd->size[1];
      b_tempd->size[0] = 1;
      b_tempd->size[1] = 1;
      emxEnsureCapacity_real_T(b_tempd, i1);
      b_tempd->data[0] = rtNaN;
    } else if (d1 - 1.0 < 1.0) {
      b_tempd->size[0] = 1;
      b_tempd->size[1] = 0;
    } else if (rtIsInf(d1 - 1.0) && (1.0 == d1 - 1.0)) {
      i1 = b_tempd->size[0] * b_tempd->size[1];
      b_tempd->size[0] = 1;
      b_tempd->size[1] = 1;
      emxEnsureCapacity_real_T(b_tempd, i1);
      b_tempd->data[0] = rtNaN;
    } else {
      i1 = b_tempd->size[0] * b_tempd->size[1];
      b_tempd->size[0] = 1;
      loop_ub = (int)floor((d1 - 1.0) - 1.0);
      b_tempd->size[1] = loop_ub + 1;
      emxEnsureCapacity_real_T(b_tempd, i1);
      for (i1 = 0; i1 <= loop_ub; i1++) {
        b_tempd->data[i1] = (double)i1 + 1.0;
      }
    }

    sparse_parenReference(Aeq1.d, Aeq1.colidx, Aeq1.rowidx, Aeq1.n, b_tempd,
                          expl_temp.d, expl_temp.colidx, expl_temp.rowidx,
                          &aeq1_temp1_m, &aeq1_temp1_n, &c_i);
    if (rtIsNaN(d1 + 1.0)) {
      i1 = b_tempd->size[0] * b_tempd->size[1];
      b_tempd->size[0] = 1;
      b_tempd->size[1] = 1;
      emxEnsureCapacity_real_T(b_tempd, i1);
      b_tempd->data[0] = rtNaN;
    } else if (Aeq1.m < appo_rep->data[(int)b_appo - 1] + 1.0) {
      b_tempd->size[0] = 1;
      b_tempd->size[1] = 0;
    } else if (rtIsInf(appo_rep->data[(int)b_appo - 1] + 1.0) && (appo_rep->
                data[(int)b_appo - 1] + 1.0 == Aeq1.m)) {
      i1 = b_tempd->size[0] * b_tempd->size[1];
      b_tempd->size[0] = 1;
      b_tempd->size[1] = 1;
      emxEnsureCapacity_real_T(b_tempd, i1);
      b_tempd->data[0] = rtNaN;
    } else if (floor(appo_rep->data[(int)b_appo - 1] + 1.0) == appo_rep->data
               [(int)b_appo - 1] + 1.0) {
      i1 = b_tempd->size[0] * b_tempd->size[1];
      b_tempd->size[0] = 1;
      loop_ub = (int)floor((double)Aeq1.m - (appo_rep->data[(int)b_appo - 1] +
        1.0));
      b_tempd->size[1] = loop_ub + 1;
      emxEnsureCapacity_real_T(b_tempd, i1);
      for (i1 = 0; i1 <= loop_ub; i1++) {
        b_tempd->data[i1] = (appo_rep->data[(int)b_appo - 1] + 1.0) + (double)i1;
      }
    } else {
      eml_float_colon(appo_rep->data[(int)b_appo - 1] + 1.0, Aeq1.m, b_tempd);
    }

    sparse_parenReference(Aeq1.d, Aeq1.colidx, Aeq1.rowidx, Aeq1.n, b_tempd,
                          aeq1_temp2.d, aeq1_temp2.colidx, aeq1_temp2.rowidx,
                          &aeq1_temp2.m, &aeq1_temp2.n, &aeq1_temp2.maxnz);
    sparse_vertcat(expl_temp.d, expl_temp.colidx, expl_temp.rowidx, aeq1_temp1_m,
                   aeq1_temp1_n, aeq1_temp2.d, aeq1_temp2.colidx,
                   aeq1_temp2.rowidx, aeq1_temp2.m, aeq1_temp2.n, &Aeq1);

    /* Aeq1=Aeq1(1:index_canc(i)-1,:); */
  }

  emxInit_real_T(&beqini, 1);
  emxInitStruct_sparse(&f_expl_temp);

  /* VVVVVVVVVVVVVVVVVVVVVVVVVVVVVV  to here   VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV */
  b_sparse_vertcat(Aeqini.d, Aeqini.colidx, Aeqini.rowidx, Aeqini.m, Aeqini.n,
                   Aeqini1.d, Aeqini1.colidx, Aeqini1.rowidx, Aeqini1.m,
                   Aeqini1.n, Aeq1.d, Aeq1.colidx, Aeq1.rowidx, Aeq1.m, Aeq1.n,
                   &f_expl_temp);

  /* Initial flows must be with a minus in front of them */
  i = beqini->size[0];
  beqini->size[0] = (int)l;
  emxEnsureCapacity_real_T(beqini, i);
  i = (int)l;
  for (b_i = 0; b_i < i; b_i++) {
    if (1 > b_i) {
      loop_ub = 0;
    } else {
      loop_ub = b_i;
    }

    delay_size[0] = loop_ub;
    if (0 <= loop_ub - 1) {
      memcpy(&delay_temp[0], &delay[0], loop_ub * sizeof(double));
    }

    b_d = b_gamma[b_i];
    beqini->data[b_i] = -(1.0 - b_d) * R0[b_i] - rt_powd_snf(1.0 - b_d,
      delay[b_i]) * V0[(int)(sum(delay_temp, delay_size) + 1.0) - 1];
  }

  i = appo_rep->size[0];
  appo_rep->size[0] = ini_P->size[0];
  emxEnsureCapacity_real_T(appo_rep, i);
  loop_ub = ini_P->size[0];
  emxFree_real_T(&ini_P);
  for (i = 0; i < loop_ub; i++) {
    appo_rep->data[i] = 0.0;
  }

  i = (int)l;
  for (b_i = 0; b_i < i; b_i++) {
    i1 = (int)(delay[b_i] + -1.0);
    if (0 <= i1 - 1) {
      if (1 > b_i) {
        loop_ub = 0;
      } else {
        loop_ub = b_i;
      }

      b_delay_size[0] = loop_ub;
      if (0 <= loop_ub - 1) {
        memcpy(&delay_data[0], &delay[0], loop_ub * sizeof(double));
      }
    }

    for (input_sizes_idx_1 = 0; input_sizes_idx_1 < i1; input_sizes_idx_1++) {
      appo_rep->data[input_sizes_idx_1] = -rt_powd_snf(1.0 - b_gamma[b_i],
        delay[b_i]) * V0[(int)(sum(delay_data, b_delay_size) + ((double)
        input_sizes_idx_1 + 2.0)) - 1];
    }
  }

  y = delay[0];
  for (k = 0; k < 6; k++) {
    y += delay[k + 1];
  }

  emxInit_real_T(&beq1, 1);
  loop_ub = (int)(varargin_1_tmp - y);
  i = beq1->size[0];
  beq1->size[0] = loop_ub;
  emxEnsureCapacity_real_T(beq1, i);
  for (i = 0; i < loop_ub; i++) {
    beq1->data[i] = 0.0;
  }

  i = b_templ->size[0];
  b_templ->size[0] = (beqini->size[0] + appo_rep->size[0]) + beq1->size[0];
  emxEnsureCapacity_real_T(b_templ, i);
  loop_ub = beqini->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_templ->data[i] = beqini->data[i];
  }

  loop_ub = appo_rep->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_templ->data[i + beqini->size[0]] = appo_rep->data[i];
  }

  loop_ub = beq1->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_templ->data[(i + beqini->size[0]) + appo_rep->size[0]] = beq1->data[i];
  }

  i = beq1->size[0];
  beq1->size[0] = b_templ->size[0];
  emxEnsureCapacity_real_T(beq1, i);
  loop_ub = b_templ->size[0];
  for (i = 0; i < loop_ub; i++) {
    beq1->data[i] = b_templ->data[i];
  }

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* A11 */
  /* H_i^n-H_i^{n-1}\le G_i^n */
  i = b_templ->size[0];
  b_templ->size[0] = (int)(N - 1.0) + 1;
  emxEnsureCapacity_real_T(b_templ, i);
  loop_ub = (int)(N - 1.0);
  for (i = 0; i < loop_ub; i++) {
    b_templ->data[i] = 1.0;
  }

  b_templ->data[(int)(N - 1.0)] = 0.0;
  c_diag(b_templ, appo);
  i = b_templ->size[0];
  b_templ->size[0] = (int)(N - 1.0);
  emxEnsureCapacity_real_T(b_templ, i);
  loop_ub = (int)(N - 1.0);
  for (i = 0; i < loop_ub; i++) {
    b_templ->data[i] = 1.0;
  }

  d_diag(b_templ, temp);
  loop_ub = appo->size[0] * appo->size[1];
  for (i = 0; i < loop_ub; i++) {
    appo->data[i] = -appo->data[i] + temp->data[i];
  }

  sparse(N * l, N * l, e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx,
         &b_i, &t2_n, &c_i);
  sparse(N * l, N * l, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  i = b_templ->size[0];
  b_templ->size[0] = (int)l;
  emxEnsureCapacity_real_T(b_templ, i);
  loop_ub = (int)l;
  for (i = 0; i < loop_ub; i++) {
    b_templ->data[i] = 1.0;
  }

  b_repmat(b_templ, N - 1.0, appo_rep);
  e_diag(appo_rep, l, temp);
  sparse(N * l, (l * l * N * M + l * N * M) + 2.0 * K * N, c_expl_temp.d,
         c_expl_temp.colidx, c_expl_temp.rowidx, &Ay_m, &Ay_n, &c_i);
  i = b_temp->size[0] * b_temp->size[1];
  b_temp->size[0] = temp->size[0];
  b_temp->size[1] = temp->size[1];
  emxEnsureCapacity_real_T(b_temp, i);
  loop_ub = temp->size[0] * temp->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_temp->data[i] = -temp->data[i];
  }

  emxInitStruct_sparse(&g_expl_temp);
  eye(l, b_r);
  c_kron(appo, b_r, r1);
  d_sparse_horzcat(b_i, t2_n, input_sizes_idx_1, b_sizes_idx_1, b_temp, Ay_m,
                   Ay_n, r1, &g_expl_temp);

  /* A12 */
  /* H_i^n-H_i^{n-1}\ge-G_i^n */
  sparse(N * l, N * l, e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx,
         &b_i, &t2_n, &c_i);
  sparse(N * l, N * l, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  i = b_templ->size[0];
  b_templ->size[0] = (int)l;
  emxEnsureCapacity_real_T(b_templ, i);
  loop_ub = (int)l;
  for (i = 0; i < loop_ub; i++) {
    b_templ->data[i] = 1.0;
  }

  b_repmat(b_templ, N - 1.0, appo_rep);
  e_diag(appo_rep, l, temp);
  sparse(N * l, (l * l * N * M + l * N * M) + 2.0 * K * N, c_expl_temp.d,
         c_expl_temp.colidx, c_expl_temp.rowidx, &Ay_m, &Ay_n, &c_i);
  eye(l, b_r);
  c_kron(appo, b_r, tempd);
  i = b_temp->size[0] * b_temp->size[1];
  b_temp->size[0] = temp->size[0];
  b_temp->size[1] = temp->size[1];
  emxEnsureCapacity_real_T(b_temp, i);
  loop_ub = temp->size[0] * temp->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_temp->data[i] = -temp->data[i];
  }

  i = b_delay_mat->size[0] * b_delay_mat->size[1];
  b_delay_mat->size[0] = tempd->size[0];
  b_delay_mat->size[1] = tempd->size[1];
  emxEnsureCapacity_real_T(b_delay_mat, i);
  loop_ub = tempd->size[0] * tempd->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_delay_mat->data[i] = -tempd->data[i];
  }

  emxInitStruct_sparse(&h_expl_temp);
  d_sparse_horzcat(b_i, t2_n, input_sizes_idx_1, b_sizes_idx_1, b_temp, Ay_m,
                   Ay_n, b_delay_mat, &h_expl_temp);

  /* VVVVVVVVVVVVVVVVVV  Changes from here   VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV */
  b_appo = (N - 1.0) * l;
  if (rtIsNaN(b_appo)) {
    i = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = 1;
    emxEnsureCapacity_real_T(b_tempd, i);
    b_tempd->data[0] = rtNaN;
  } else if (b_appo < 1.0) {
    b_tempd->size[0] = 1;
    b_tempd->size[1] = 0;
  } else if (rtIsInf(b_appo) && (1.0 == b_appo)) {
    i = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = 1;
    emxEnsureCapacity_real_T(b_tempd, i);
    b_tempd->data[0] = rtNaN;
  } else {
    i = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    loop_ub = (int)floor(b_appo - 1.0);
    b_tempd->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(b_tempd, i);
    for (i = 0; i <= loop_ub; i++) {
      b_tempd->data[i] = (double)i + 1.0;
    }
  }

  sparse_parenReference(g_expl_temp.d, g_expl_temp.colidx, g_expl_temp.rowidx,
                        g_expl_temp.n, b_tempd, aeq1_temp2.d, aeq1_temp2.colidx,
                        aeq1_temp2.rowidx, &aeq1_temp2.m, &aeq1_temp2.n,
                        &aeq1_temp2.maxnz);
  if (rtIsNaN(b_appo)) {
    i = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = 1;
    emxEnsureCapacity_real_T(b_tempd, i);
    b_tempd->data[0] = rtNaN;
  } else if (b_appo < 1.0) {
    b_tempd->size[0] = 1;
    b_tempd->size[1] = 0;
  } else if (rtIsInf(b_appo) && (1.0 == b_appo)) {
    i = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = 1;
    emxEnsureCapacity_real_T(b_tempd, i);
    b_tempd->data[0] = rtNaN;
  } else {
    i = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    loop_ub = (int)floor(b_appo - 1.0);
    b_tempd->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(b_tempd, i);
    for (i = 0; i <= loop_ub; i++) {
      b_tempd->data[i] = (double)i + 1.0;
    }
  }

  sparse_parenReference(h_expl_temp.d, h_expl_temp.colidx, h_expl_temp.rowidx,
                        h_expl_temp.n, b_tempd, Aeq1.d, Aeq1.colidx, Aeq1.rowidx,
                        &Aeq1.m, &Aeq1.n, &Aeq1.maxnz);

  /* VVVVVVVVVVVVVVVVVVVVVVVVVVVVVV  to here   VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV */
  i = appo_rep->size[0];
  appo_rep->size[0] = (int)b_appo;
  emxEnsureCapacity_real_T(appo_rep, i);
  loop_ub_tmp = (int)(l * (N - 1.0));
  for (i = 0; i < loop_ub_tmp; i++) {
    appo_rep->data[i] = 1.0;
  }

  i = beqini->size[0];
  beqini->size[0] = loop_ub_tmp;
  emxEnsureCapacity_real_T(beqini, i);
  for (i = 0; i < loop_ub_tmp; i++) {
    beqini->data[i] = -1.0;
  }

  /* Aini1 */
  /* initial gate opening conditions: The gate opening ratio cannot change at the beginning if the gate is not operated  */
  sparse(l, l * N, e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx, &b_i,
         &t2_n, &c_i);
  sparse(l, l * N, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  eye(l, temp);
  sparse(l, (N - 1.0) * l, c_expl_temp.d, c_expl_temp.colidx, c_expl_temp.rowidx,
         &Ay_m, &Ay_n, &c_i);
  sparse(l, (l * l * N * M + l * N * M) + 2.0 * K * N, b_expl_temp.d,
         b_expl_temp.colidx, b_expl_temp.rowidx, &k, &sizes_idx_1, &c_i);
  sparse(l, (N - 1.0) * l, A102.d, A102.colidx, A102.rowidx, &A102.m, &A102.n,
         &A102.maxnz);
  i = b_temp->size[0] * b_temp->size[1];
  b_temp->size[0] = temp->size[0];
  b_temp->size[1] = temp->size[1];
  emxEnsureCapacity_real_T(b_temp, i);
  loop_ub = temp->size[0] * temp->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_temp->data[i] = -temp->data[i];
  }

  eye(l, b_r);
  e_sparse_horzcat(b_i, t2_n, input_sizes_idx_1, b_sizes_idx_1, b_temp, Ay_m,
                   Ay_n, k, sizes_idx_1, b_r, A102.m, A102.n, &expl_temp);
  sparse(l, l * N, e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx, &b_i,
         &t2_n, &c_i);
  sparse(l, l * N, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  eye(l, temp);
  sparse(l, (N - 1.0) * l, c_expl_temp.d, c_expl_temp.colidx, c_expl_temp.rowidx,
         &Ay_m, &Ay_n, &c_i);
  sparse(l, (l * l * N * M + l * N * M) + 2.0 * K * N, b_expl_temp.d,
         b_expl_temp.colidx, b_expl_temp.rowidx, &k, &sizes_idx_1, &c_i);
  eye(l, b_r);
  sparse(l, (N - 1.0) * l, A102.d, A102.colidx, A102.rowidx, &A102.m, &A102.n,
         &A102.maxnz);
  i = b_temp->size[0] * b_temp->size[1];
  b_temp->size[0] = temp->size[0];
  b_temp->size[1] = temp->size[1];
  emxEnsureCapacity_real_T(b_temp, i);
  loop_ub = temp->size[0] * temp->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_temp->data[i] = -temp->data[i];
  }

  i = b_delay_mat->size[0] * b_delay_mat->size[1];
  b_delay_mat->size[0] = b_r->size[0];
  b_delay_mat->size[1] = b_r->size[1];
  emxEnsureCapacity_real_T(b_delay_mat, i);
  loop_ub = b_r->size[0] * b_r->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_delay_mat->data[i] = -b_r->data[i];
  }

  e_sparse_horzcat(b_i, t2_n, input_sizes_idx_1, b_sizes_idx_1, b_temp, Ay_m,
                   Ay_n, k, sizes_idx_1, b_delay_mat, A102.m, A102.n, &Aeqini);
  sparse_vertcat(expl_temp.d, expl_temp.colidx, expl_temp.rowidx, expl_temp.m,
                 expl_temp.n, aeq1_temp2.d, aeq1_temp2.colidx, aeq1_temp2.rowidx,
                 aeq1_temp2.m, aeq1_temp2.n, &c_expl_temp);
  i = g_expl_temp.d->size[0];
  g_expl_temp.d->size[0] = c_expl_temp.d->size[0];
  emxEnsureCapacity_real_T(g_expl_temp.d, i);
  loop_ub = c_expl_temp.d->size[0];
  for (i = 0; i < loop_ub; i++) {
    g_expl_temp.d->data[i] = c_expl_temp.d->data[i];
  }

  i = g_expl_temp.colidx->size[0];
  g_expl_temp.colidx->size[0] = c_expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(g_expl_temp.colidx, i);
  loop_ub = c_expl_temp.colidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    g_expl_temp.colidx->data[i] = c_expl_temp.colidx->data[i];
  }

  i = g_expl_temp.rowidx->size[0];
  g_expl_temp.rowidx->size[0] = c_expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(g_expl_temp.rowidx, i);
  loop_ub = c_expl_temp.rowidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    g_expl_temp.rowidx->data[i] = c_expl_temp.rowidx->data[i];
  }

  A11_m = c_expl_temp.m;
  A11_n = c_expl_temp.n;
  sparse_vertcat(Aeqini.d, Aeqini.colidx, Aeqini.rowidx, Aeqini.m, Aeqini.n,
                 Aeq1.d, Aeq1.colidx, Aeq1.rowidx, Aeq1.m, Aeq1.n, &c_expl_temp);
  i = h_expl_temp.d->size[0];
  h_expl_temp.d->size[0] = c_expl_temp.d->size[0];
  emxEnsureCapacity_real_T(h_expl_temp.d, i);
  loop_ub = c_expl_temp.d->size[0];
  for (i = 0; i < loop_ub; i++) {
    h_expl_temp.d->data[i] = c_expl_temp.d->data[i];
  }

  i = h_expl_temp.colidx->size[0];
  h_expl_temp.colidx->size[0] = c_expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(h_expl_temp.colidx, i);
  loop_ub = c_expl_temp.colidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    h_expl_temp.colidx->data[i] = c_expl_temp.colidx->data[i];
  }

  i = h_expl_temp.rowidx->size[0];
  h_expl_temp.rowidx->size[0] = c_expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(h_expl_temp.rowidx, i);
  loop_ub = c_expl_temp.rowidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    h_expl_temp.rowidx->data[i] = c_expl_temp.rowidx->data[i];
  }

  A12_m = c_expl_temp.m;
  A12_n = c_expl_temp.n;

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* A2 */
  /* V_i^n\le  min\{c_i,r_n \}  */
  i = b_templ->size[0];
  b_templ->size[0] = (int)varargin_1_tmp;
  emxEnsureCapacity_real_T(b_templ, i);
  for (i = 0; i < index_canc_idx_0; i++) {
    b_templ->data[i] = 1.0;
  }

  spdiags(b_templ, N * l, N * l, e_expl_temp.d, e_expl_temp.colidx,
          e_expl_temp.rowidx, &b_i, &t2_n);
  sparse(N * l, ((2.0 * N * l + l * l * N * M) + l * N * M) + 2.0 * K * N,
         d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  b_diag(l,c, psi_temp);
  eye(N, r1);
  b_kron(l,r1, psi_temp, b_r);
  i = b_delay_mat->size[0] * b_delay_mat->size[1];
  b_delay_mat->size[0] = b_r->size[0];
  b_delay_mat->size[1] = b_r->size[1];
  emxEnsureCapacity_real_T(b_delay_mat, i);
  loop_ub = b_r->size[0] * b_r->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_delay_mat->data[i] = -b_r->data[i];
  }

  emxInit_int8_T(&b21, 1);
  emxInitStruct_sparse(&i_expl_temp);
  f_sparse_horzcat(e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx, b_i,
                   t2_n, input_sizes_idx_1, b_sizes_idx_1, b_delay_mat,
                   &i_expl_temp);
  i = b21->size[0];
  b21->size[0] = (int)varargin_1_tmp;
  emxEnsureCapacity_int8_T(b21, i);
  for (i = 0; i < index_canc_idx_0; i++) {
    b21->data[i] = 0;
  }

  i = b_templ->size[0];
  b_templ->size[0] = (int)varargin_1_tmp;
  emxEnsureCapacity_real_T(b_templ, i);
  for (i = 0; i < index_canc_idx_0; i++) {
    b_templ->data[i] = 1.0;
  }

  spdiags(b_templ, N * l, N * l, expl_temp.d, expl_temp.colidx, expl_temp.rowidx,
          &aeq1_temp1_m, &aeq1_temp1_n);
  loop_ub = expl_temp.d->size[0];
  for (i = 0; i < loop_ub; i++) {
    expl_temp.d->data[i] = -expl_temp.d->data[i];
  }

  sparse(N * l, ((2.0 * N * l + l * l * N * M) + l * N * M) + 2.0 * K * N,
         e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx, &b_i, &t2_n,
         &c_i);
  for (i = 0; i < l; i++) {
    delay_temp[i] = c[i] * rho[i];
  }

  b_diag(l,delay_temp, psi_temp);
  emxInitStruct_sparse(&A22);
  emxInit_int8_T(&b22, 1);
  eye(N, b_r);
  b_kron(l,b_r, psi_temp, r1);
  f_sparse_horzcat(expl_temp.d, expl_temp.colidx, expl_temp.rowidx, aeq1_temp1_m,
                   aeq1_temp1_n, b_i, t2_n, r1, &A22);
  i = b22->size[0];
  b22->size[0] = (int)varargin_1_tmp;
  emxEnsureCapacity_int8_T(b22, i);
  for (i = 0; i < index_canc_idx_0; i++) {
    b22->data[i] = 0;
  }

  i = b_tempd->size[0] * b_tempd->size[1];
  b_tempd->size[0] = 1;
  b_tempd->size[1] = (int)(l - 1.0) + 1;
  emxEnsureCapacity_real_T(b_tempd, i);
  b_tempd->data[0] = 1.0;
  loop_ub = (int)(l - 1.0);
  for (i = 0; i < loop_ub; i++) {
    b_tempd->data[i + 1] = 0.0;
  }

  eye(N, r1);
  kron(r1, b_tempd, b_r);
  sparse(N, ((2.0 * l * N + l * l * N * M) + l * N * M) + 2.0 * K * N,
         e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx, &b_i, &t2_n,
         &c_i);
  i = b_tempd->size[0] * b_tempd->size[1];
  b_tempd->size[0] = 1;
  b_tempd->size[1] = (int)(l - 1.0) + 1;
  emxEnsureCapacity_real_T(b_tempd, i);
  b_tempd->data[0] = c[0];
  loop_ub = (int)(l - 1.0);
  for (i = 0; i < loop_ub; i++) {
    b_tempd->data[i + 1] = 0.0;
  }

  eye(N, r1);
  kron(r1, b_tempd, tempd);
  i = b_delay_mat->size[0] * b_delay_mat->size[1];
  b_delay_mat->size[0] = tempd->size[0];
  b_delay_mat->size[1] = tempd->size[1];
  emxEnsureCapacity_real_T(b_delay_mat, i);
  loop_ub = tempd->size[0] * tempd->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_delay_mat->data[i] = -tempd->data[i];
  }

  emxInitStruct_sparse(&A2eqini);
  g_sparse_horzcat(b_r, b_i, t2_n, b_delay_mat, &A2eqini);

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* Aeq2 */
  /* \sum_{i=1}^l\ G_i^n= 0, for n \in t^p */
  i = Aeq2->size[0] * Aeq2->size[1];
  Aeq2->size[0] = (int)tp_size;
  b_d = varargin_1_tmp * M;
  d1 = 2.0 * K * N;
  y = rt_powd_snf(l, 2.0) * N * M;
  a = ((3.0 * l * N + y) + b_d) + d1;
  i1 = (int)(a + varargin_1_tmp);
  Aeq2->size[1] = i1;
  emxEnsureCapacity_real_T(Aeq2, i);
  loop_ub = (int)tp_size * i1;
  for (i = 0; i < loop_ub; i++) {
    Aeq2->data[i] = 0.0;
  }

  /*  for G */
  i = (int)tp_size;
  for (b_i = 0; b_i < i; b_i++) {
    d2 = 2.0 * l * N;
    b_appo = tp[b_i];
    b_y = (d2 + (b_appo - 1.0) * l) + 1.0;
    d2 += b_appo * l;
    if (b_y > d2) {
      Ay_n = 0;
      Ay_m = 0;
    } else {
      Ay_n = (int)b_y - 1;
      Ay_m = (int)d2;
    }

    loop_ub = Ay_m - Ay_n;
    for (Ay_m = 0; Ay_m < loop_ub; Ay_m++) {
      Aeq2->data[b_i + Aeq2->size[0] * (Ay_n + Ay_m)] = 1.0;
    }
  }

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* A3 */
  /*  \sum_{i=1}^{l}\sum_{n=1}^N E_i^{n,m}\le 1,\qquad \forall m=1,\dots,M   */
  sparse(M, 3.0 * l * N + N * M * l * l, e_expl_temp.d, e_expl_temp.colidx,
         e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  i = b_tempd->size[0] * b_tempd->size[1];
  b_tempd->size[0] = 1;
  b_tempd->size[1] = (int)l;
  emxEnsureCapacity_real_T(b_tempd, i);
  loop_ub = (int)l;
  for (i = 0; i < loop_ub; i++) {
    b_tempd->data[i] = 1.0;
  }

  emxInitStruct_sparse(&A3);
  sparse(M, 2.0 * K * N, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  sparse(M, l * N, c_expl_temp.d, c_expl_temp.colidx, c_expl_temp.rowidx, &Ay_m,
         &Ay_n, &c_i);
  eye(M, b_r);
  kron(b_r, b_tempd, r1);
  c_repmat(r1, N, b_r);
  h_sparse_horzcat(b_i, t2_n, b_r, input_sizes_idx_1, b_sizes_idx_1, Ay_m, Ay_n,
                   &A3);

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* Aeq3 */
  /* G_i^n- \sum_{m=1}^M E_i^{n,m}&=0 */
  sparse(l * N, 2.0 * l * N, e_expl_temp.d, e_expl_temp.colidx,
         e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  sparse(l * N, rt_powd_snf(l, 2.0) * N * M, d_expl_temp.d, d_expl_temp.colidx,
         d_expl_temp.rowidx, &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  eye(N, b_r);
  eye(l, r1);
  c_repmat(r1, M, b_delay_mat);
  c_kron(b_r, b_delay_mat, tempd);
  sparse(l * N, 2.0 * K * N, c_expl_temp.d, c_expl_temp.colidx,
         c_expl_temp.rowidx, &Ay_m, &Ay_n, &c_i);
  sparse(l * N, l * N, b_expl_temp.d, b_expl_temp.colidx, b_expl_temp.rowidx, &k,
         &sizes_idx_1, &c_i);
  i = b_delay_mat->size[0] * b_delay_mat->size[1];
  b_delay_mat->size[0] = tempd->size[0];
  b_delay_mat->size[1] = tempd->size[1];
  emxEnsureCapacity_real_T(b_delay_mat, i);
  loop_ub = tempd->size[0] * tempd->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_delay_mat->data[i] = -tempd->data[i];
  }

  emxInit_int8_T(&beq3, 1);
  emxInitStruct_sparse(&j_expl_temp);
  eye(varargin_1_tmp, b_r);
  i_sparse_horzcat(b_i, t2_n, b_r, input_sizes_idx_1, b_sizes_idx_1, b_delay_mat,
                   Ay_m, Ay_n, k, sizes_idx_1, &j_expl_temp);
  i = beq3->size[0];
  beq3->size[0] = (int)varargin_1_tmp;
  emxEnsureCapacity_int8_T(beq3, i);
  for (i = 0; i < index_canc_idx_0; i++) {
    beq3->data[i] = 0;
  }

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* A4 */
  /* E_i^{m-1,n}-\sum_{j=1}^l\sum_{p=n}^N F_{ij}^{m,p}\le0 */
  index_canc_idx_0_tmp = (M - 1.0) * l;
  k = (int)index_canc_idx_0_tmp;
  i = b_templ->size[0];
  b_templ->size[0] = (int)index_canc_idx_0_tmp;
  emxEnsureCapacity_real_T(b_templ, i);
  for (i = 0; i < k; i++) {
    b_templ->data[i] = 1.0;
  }

  e_diag(b_templ, l, appo);
  b_varargin_1 = M * l;
  i = appo_tmp->size[0] * appo_tmp->size[1];
  appo_tmp->size[0] = 1;
  loop_ub = (int)floor(b_varargin_1 - (index_canc_idx_0_tmp + 1.0));
  appo_tmp->size[1] = loop_ub + 1;
  emxEnsureCapacity_int32_T(appo_tmp, i);
  for (i = 0; i <= loop_ub; i++) {
    appo_tmp->data[i] = (int)((index_canc_idx_0_tmp + 1.0) + (double)i);
  }

  c_nullAssignment(appo, appo_tmp);
  i = temp->size[0] * temp->size[1];
  temp->size[0] = (int)N;
  temp->size[1] = (int)N;
  emxEnsureCapacity_real_T(temp, i);
  loop_ub = (int)N * (int)N;
  for (i = 0; i < loop_ub; i++) {
    temp->data[i] = 1.0;
  }

  triu(temp);
  i = b_tempd->size[0] * b_tempd->size[1];
  b_tempd->size[0] = 1;
  b_tempd->size[1] = (int)l;
  emxEnsureCapacity_real_T(b_tempd, i);
  loop_ub = (int)l;
  for (i = 0; i < loop_ub; i++) {
    b_tempd->data[i] = 1.0;
  }

  c_kron(temp, appo, b_r);
  kron(b_r, b_tempd, appo);
  loop_ub = appo->size[0] * appo->size[1];
  for (i = 0; i < loop_ub; i++) {
    appo->data[i] = -appo->data[i];
  }

  emxInitStruct_sparse(&A4);
  emxInit_int8_T(&b4, 1);
  sparse(l * N * (M - 1.0), 3.0 * l * N, e_expl_temp.d, e_expl_temp.colidx,
         e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  eye(b_varargin_1, b_r);
  c_nullAssignment(b_r, appo_tmp);
  sparse(l * N * (M - 1.0), 2.0 * K * N, d_expl_temp.d, d_expl_temp.colidx,
         d_expl_temp.rowidx, &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  sparse(l * N * (M - 1.0), l * N, c_expl_temp.d, c_expl_temp.colidx,
         c_expl_temp.rowidx, &Ay_m, &Ay_n, &c_i);
  eye(N, r1);
  c_kron(r1, b_r, b_delay_mat);
  j_sparse_horzcat(b_i, t2_n, appo, b_delay_mat, input_sizes_idx_1,
                   b_sizes_idx_1, Ay_m, Ay_n, &A4);
  loop_ub = (int)(varargin_1_tmp * (M - 1.0));
  i = b4->size[0];
  b4->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(b4, i);
  for (i = 0; i < loop_ub; i++) {
    b4->data[i] = 0;
  }

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* Aeq4 */
  /* E_j^{m,n}= \sum_{i=1}^l F_{ij}^{m,n} */
  sparse(l * N * M, 3.0 * l * N, e_expl_temp.d, e_expl_temp.colidx,
         e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  d2 = M * N;
  eye(d2, r1);
  eye(l, b_r);
  c_repmat(b_r, l, b_delay_mat);
  c_kron(r1, b_delay_mat, b_r);
  index_canc_idx_0_tmp = b_varargin_1 * N;
  k = (int)index_canc_idx_0_tmp;
  i = b_templ->size[0];
  b_templ->size[0] = (int)index_canc_idx_0_tmp;
  emxEnsureCapacity_real_T(b_templ, i);
  for (i = 0; i < k; i++) {
    b_templ->data[i] = 1.0;
  }

  spdiags(b_templ, l * M * N, l * M * N, d_expl_temp.d, d_expl_temp.colidx,
          d_expl_temp.rowidx, &input_sizes_idx_1, &b_sizes_idx_1);
  sparse(l * N * M, 2.0 * K * N, c_expl_temp.d, c_expl_temp.colidx,
         c_expl_temp.rowidx, &Ay_m, &Ay_n, &c_i);
  sparse(l * N * M, l * N, b_expl_temp.d, b_expl_temp.colidx, b_expl_temp.rowidx,
         &k, &sizes_idx_1, &c_i);
  i = b_delay_mat->size[0] * b_delay_mat->size[1];
  b_delay_mat->size[0] = b_r->size[0];
  b_delay_mat->size[1] = b_r->size[1];
  emxEnsureCapacity_real_T(b_delay_mat, i);
  loop_ub = b_r->size[0] * b_r->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_delay_mat->data[i] = -b_r->data[i];
  }

  emxInit_int8_T(&beq4, 1);
  emxInitStruct_sparse(&k_expl_temp);
  k_sparse_horzcat(b_i, t2_n, b_delay_mat, d_expl_temp.d, d_expl_temp.colidx,
                   d_expl_temp.rowidx, input_sizes_idx_1, b_sizes_idx_1, Ay_m,
                   Ay_n, k, sizes_idx_1, &k_expl_temp);
  i = beq4->size[0];
  beq4->size[0] = (int)b_d;
  emxEnsureCapacity_int8_T(beq4, i);
  loop_ub_tmp = (int)(l * N * M);
  for (i = 0; i < loop_ub_tmp; i++) {
    beq4->data[i] = 0;
  }

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* A5 */
  /* \sum_{j=1}^l\sum_{n=1}^N nE_j^{n,m}-\sum_{j=1}^l\sum_{n=1}^N nE_j^{n,m-1}\ge  \sum_{i=1}^l\sum_{j=1}^l\sum_{n=1}^N \psi_{ij}F_{ij}^{n,m} */
  for (b_i = 0; b_i < l*l; b_i++) {
    b_appo = psi[b_i];
    psi_temp[b_i] = b_appo;
    if (b_appo <= 1.0) {
      psi_temp[b_i] = 0.0;
    }
  }

  if (N < 1.0) {
    b_tempd->size[0] = 1;
    b_tempd->size[1] = 0;
  } else if (rtIsInf(N) && (1.0 == N)) {
    i = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = 1;
    emxEnsureCapacity_real_T(b_tempd, i);
    b_tempd->data[0] = rtNaN;
  } else {
    i = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    loop_ub = (int)floor(N - 1.0);
    b_tempd->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(b_tempd, i);
    for (i = 0; i <= loop_ub; i++) {
      b_tempd->data[i] = (double)i + 1.0;
    }
  }

  sparse(M, 3.0 * l * N, e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx,
         &b_i, &t2_n, &c_i);
  i = b_templ->size[0];
  b_templ->size[0] = (int)(M - 1.0);
  emxEnsureCapacity_real_T(b_templ, i);
  loop_ub = (int)(M - 1.0);
  for (i = 0; i < loop_ub; i++) {
    b_templ->data[i] = 1.0;
  }

  d_diag(b_templ, b_r);
  c_repmat(b_r, N, r1);
  d_kron(l,r1, psi_temp, temp);
  i = b_templ->size[0];
  b_templ->size[0] = (int)M;
  emxEnsureCapacity_real_T(b_templ, i);
  loop_ub = (int)M;
  for (i = 0; i < loop_ub; i++) {
    b_templ->data[i] = 1.0;
  }

  f_diag(b_templ, b_r);
  i = b_templ->size[0];
  b_templ->size[0] = (int)(M - 1.0);
  emxEnsureCapacity_real_T(b_templ, i);
  loop_ub = (int)(M - 1.0);
  for (i = 0; i < loop_ub; i++) {
    b_templ->data[i] = 1.0;
  }

  d_diag(b_templ, tempd);
  i = b_delay_mat->size[0] * b_delay_mat->size[1];
  b_delay_mat->size[0] = b_r->size[0];
  b_delay_mat->size[1] = b_r->size[1];
  emxEnsureCapacity_real_T(b_delay_mat, i);
  loop_ub = b_r->size[0] * b_r->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_delay_mat->data[i] = b_r->data[i] - tempd->data[i];
  }

  emxFree_real_T(&tempd);
  e_kron(b_tempd, b_delay_mat, b_r);
  i = b_tempd->size[0] * b_tempd->size[1];
  b_tempd->size[0] = 1;
  b_tempd->size[1] = (int)l;
  emxEnsureCapacity_real_T(b_tempd, i);
  loop_ub = (int)l;
  for (i = 0; i < loop_ub; i++) {
    b_tempd->data[i] = 1.0;
  }

  emxInitStruct_sparse(&A5);
  emxInit_real_T(&t4_d, 1);
  emxInit_int32_T(&t4_colidx, 1);
  emxInit_int32_T(&t4_rowidx, 1);
  sparse(M, 2.0 * K * N, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  sparse(M, l * N, c_expl_temp.d, c_expl_temp.colidx, c_expl_temp.rowidx, &Ay_m,
         &Ay_n, &c_i);
  kron(b_r, b_tempd, r1);
  j_sparse_horzcat(b_i, t2_n, temp, r1, input_sizes_idx_1, b_sizes_idx_1, Ay_m,
                   Ay_n, &A5);
  c_sparse(((4.0 * l * N + rt_powd_snf(l, 2.0) * N * M) + l * N * M) + 2.0 * K *
           N, t4_d, t4_colidx, t4_rowidx, &input_sizes_idx_1);
  sparse_parenAssign(&A5, input_sizes_idx_1, M);

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* A6 */
  /*  \sum_{i=1}^l\sum_{j=1}^l\sum_{m=1}^M \psi_{ij}F_{ij}^{n,m}\le 1 */
  for (b_i = 0; b_i < l*l; b_i++) {
    b_appo = psi[b_i];
    psi_temp[b_i] = b_appo;
    if (b_appo > 1.0) {
      psi_temp[b_i] = 0.0;
    }
  }

  emxInitStruct_sparse(&A6);
  sparse(N, 3.0 * l * N, e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx,
         &b_i, &t2_n, &c_i);
  sparse(N, N * M * l, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  sparse(N, 2.0 * K * N, c_expl_temp.d, c_expl_temp.colidx, c_expl_temp.rowidx,
         &Ay_m, &Ay_n, &c_i);
  sparse(N, l * N, b_expl_temp.d, b_expl_temp.colidx, b_expl_temp.rowidx, &k,
         &sizes_idx_1, &c_i);
  eye(N, b_r);
  d_repmat(l,psi_temp, M, b_tempd);
  kron(b_r, b_tempd, r1);
  k_sparse_horzcat(b_i, t2_n, r1, d_expl_temp.d, d_expl_temp.colidx,
                   d_expl_temp.rowidx, input_sizes_idx_1, b_sizes_idx_1, Ay_m,
                   Ay_n, k, sizes_idx_1, &A6);

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* A7 */
  /*  \sum_{n=1}^N S_k^n=1 */
  sparse(K, (3.0 * l * N + N * M * l) + N * M * l * l, e_expl_temp.d,
         e_expl_temp.colidx, e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  i = b_templ->size[0];
  b_templ->size[0] = (int)K;
  emxEnsureCapacity_real_T(b_templ, i);
  loop_ub = (int)K;
  for (i = 0; i < loop_ub; i++) {
    b_templ->data[i] = 1.0;
  }

  emxInitStruct_sparse(&A7);
  spdiags(b_templ, K, K, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
          &input_sizes_idx_1, &b_sizes_idx_1);
  sparse_repmat(d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
                input_sizes_idx_1, b_sizes_idx_1, N, &c_expl_temp);
  sparse(K, K * N, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  sparse(K, l * N, b_expl_temp.d, b_expl_temp.colidx, b_expl_temp.rowidx, &k,
         &sizes_idx_1, &c_i);
  l_sparse_horzcat(b_i, t2_n, c_expl_temp.d, c_expl_temp.colidx,
                   c_expl_temp.rowidx, c_expl_temp.m, c_expl_temp.n,
                   d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
                   input_sizes_idx_1, b_sizes_idx_1, k, sizes_idx_1, &A7);

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* Aeq5 */
  /*  S_k^1= D_k^1 */
  sparse(K, K * (N - 1.0), e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx,
         &b_i, &t2_n, &c_i);
  eye(K, b_r);
  m_sparse_horzcat(b_r, b_i, t2_n, &aeq1_temp2);
  i = expl_temp.d->size[0];
  expl_temp.d->size[0] = aeq1_temp2.d->size[0];
  emxEnsureCapacity_real_T(expl_temp.d, i);
  loop_ub = aeq1_temp2.d->size[0];
  for (i = 0; i < loop_ub; i++) {
    expl_temp.d->data[i] = aeq1_temp2.d->data[i];
  }

  i = expl_temp.colidx->size[0];
  expl_temp.colidx->size[0] = aeq1_temp2.colidx->size[0];
  emxEnsureCapacity_int32_T(expl_temp.colidx, i);
  loop_ub = aeq1_temp2.colidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    expl_temp.colidx->data[i] = aeq1_temp2.colidx->data[i];
  }

  i = expl_temp.rowidx->size[0];
  expl_temp.rowidx->size[0] = aeq1_temp2.rowidx->size[0];
  emxEnsureCapacity_int32_T(expl_temp.rowidx, i);
  loop_ub = aeq1_temp2.rowidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    expl_temp.rowidx->data[i] = aeq1_temp2.rowidx->data[i];
  }

  loop_ub = expl_temp.d->size[0];
  for (i = 0; i < loop_ub; i++) {
    expl_temp.d->data[i] = -expl_temp.d->data[i];
  }

  emxInitStruct_sparse(&Aeq5);
  sparse(K, (3.0 * l * N + N * M * l) + N * M * l * l, e_expl_temp.d,
         e_expl_temp.colidx, e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  sparse(K, l * N, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  l_sparse_horzcat(b_i, t2_n, aeq1_temp2.d, aeq1_temp2.colidx, aeq1_temp2.rowidx,
                   aeq1_temp2.m, aeq1_temp2.n, expl_temp.d, expl_temp.colidx,
                   expl_temp.rowidx, aeq1_temp2.m, aeq1_temp2.n,
                   input_sizes_idx_1, b_sizes_idx_1, &Aeq5);

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* A8 */
  /*  S_k^n=0\rightarrow D_k^{n-1}-D_k^n\ge0 */
  b_appo = K * (N - 1.0);
  k = (int)b_appo;
  i = b_templ->size[0];
  b_templ->size[0] = (int)b_appo;
  emxEnsureCapacity_real_T(b_templ, i);
  for (i = 0; i < k; i++) {
    b_templ->data[i] = 1.0;
  }

  e_diag(b_templ, K, appo);
  loop_ub = appo->size[0] * appo->size[1];
  for (i = 0; i < loop_ub; i++) {
    appo->data[i] = -appo->data[i];
  }

  i = appo_tmp->size[0] * appo_tmp->size[1];
  appo_tmp->size[0] = 1;
  loop_ub = (int)floor(varargin_1 - (b_appo + 1.0));
  appo_tmp->size[1] = loop_ub + 1;
  emxEnsureCapacity_int32_T(appo_tmp, i);
  for (i = 0; i <= loop_ub; i++) {
    appo_tmp->data[i] = (int)((b_appo + 1.0) + (double)i);
  }

  c_nullAssignment(appo, appo_tmp);
  k = (int)(K * (N - 1.0));
  i = b_templ->size[0];
  b_templ->size[0] = k;
  emxEnsureCapacity_real_T(b_templ, i);
  for (i = 0; i < k; i++) {
    b_templ->data[i] = 1.0;
  }

  e_diag(b_templ, -K, temp);
  sparse(K * (N - 1.0), (3.0 * l * N + N * M * l) + N * M * l * l, e_expl_temp.d,
         e_expl_temp.colidx, e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  eye(varargin_1, b_r);
  loop_ub = b_r->size[0] * b_r->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_r->data[i] -= temp->data[i];
  }

  i = appo_tmp->size[0] * appo_tmp->size[1];
  appo_tmp->size[0] = 1;
  loop_ub = (int)floor(K - 1.0);
  appo_tmp->size[1] = loop_ub + 1;
  emxEnsureCapacity_int32_T(appo_tmp, i);
  for (i = 0; i <= loop_ub; i++) {
    appo_tmp->data[i] = i + 1;
  }

  emxInit_int8_T(&b8, 1);
  c_nullAssignment(b_r, appo_tmp);
  sparse(K * (N - 1.0), l * N, d_expl_temp.d, d_expl_temp.colidx,
         d_expl_temp.rowidx, &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  n_sparse_horzcat(b_i, t2_n, appo, b_r, input_sizes_idx_1, b_sizes_idx_1, &Aeq1);

  /* A82=[sparse(K*(N-1),3*l*N+N*M*l+N*M*l*l),-appo,-appo1,sparse(K*(N-1),l*N)]; */
  i = b8->size[0];
  b8->size[0] = k;
  emxEnsureCapacity_int8_T(b8, i);
  emxFree_int32_T(&appo_tmp);
  emxFree_real_T(&appo);
  for (i = 0; i < k; i++) {
    b8->data[i] = 0;
  }

  /* b82=sparse(K*(N-1),1); */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* A9 */
  /*  \sum_{n=1}^N nS_k^n+\sum_{n=1}^N D_k^n&\le N */
  if (N < 1.0) {
    b_tempd->size[0] = 1;
    b_tempd->size[1] = 0;
  } else if (rtIsInf(N) && (1.0 == N)) {
    i = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = 1;
    emxEnsureCapacity_real_T(b_tempd, i);
    b_tempd->data[0] = rtNaN;
  } else {
    i = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = (int)floor(N - 1.0) + 1;
    emxEnsureCapacity_real_T(b_tempd, i);
    loop_ub = (int)floor(N - 1.0);
    for (i = 0; i <= loop_ub; i++) {
      b_tempd->data[i] = (double)i + 1.0;
    }
  }

  sparse(K, (3.0 * l * N + N * M * l) + N * M * l * l, e_expl_temp.d,
         e_expl_temp.colidx, e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  eye(K, r1);
  e_kron(b_tempd, r1, b_r);
  i = b_tempd->size[0] * b_tempd->size[1];
  b_tempd->size[0] = 1;
  b_tempd->size[1] = (int)N;
  emxEnsureCapacity_real_T(b_tempd, i);
  loop_ub = (int)N;
  for (i = 0; i < loop_ub; i++) {
    b_tempd->data[i] = 1.0;
  }

  sparse(K, l * N, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  eye(K, r1);
  e_kron(b_tempd, r1, b_delay_mat);
  n_sparse_horzcat(b_i, t2_n, b_r, b_delay_mat, input_sizes_idx_1, b_sizes_idx_1,
                   &Aeqini);

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* A10 */
  /* \epsilon_k q_k d_k\le\sum_{n=1}^N \Delta TD_k^n\le q_k d_k\qquad \forall k=1,\dots,K,  */
  sparse(K, ((3.0 * l * N + N * M * l) + N * M * l * l) + K * N, e_expl_temp.d,
         e_expl_temp.colidx, e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  i = b_tempd->size[0] * b_tempd->size[1];
  b_tempd->size[0] = 1;
  b_tempd->size[1] = (int)N;
  emxEnsureCapacity_real_T(b_tempd, i);
  loop_ub = (int)N;
  for (i = 0; i < loop_ub; i++) {
    b_tempd->data[i] = 1.0;
  }

  emxInit_real_T(&r2, 2);
  sparse(K, l * N, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  g_diag(K,q, dv);
  f_kron(K,b_tempd, dv, r2);
  o_sparse_horzcat(K,b_i, t2_n, r2, input_sizes_idx_1, b_sizes_idx_1, &Aeqini1);
  sparse(K, (3.0 * l * N + N * M * l) + N * M * l * l, e_expl_temp.d,
         e_expl_temp.colidx, e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  for (i = 0; i < K; i++) {
    b_eps[i] = -eps[i] * d[i] * q[i];
  }

  b_spdiags(K,b_eps, K, K, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
            &input_sizes_idx_1, &b_sizes_idx_1);
  sparse_repmat(d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
                input_sizes_idx_1, b_sizes_idx_1, N, &e_expl_temp);
  i = c_expl_temp.d->size[0];
  c_expl_temp.d->size[0] = e_expl_temp.d->size[0];
  emxEnsureCapacity_real_T(c_expl_temp.d, i);
  loop_ub = e_expl_temp.d->size[0];
  for (i = 0; i < loop_ub; i++) {
    c_expl_temp.d->data[i] = e_expl_temp.d->data[i];
  }

  i = c_expl_temp.colidx->size[0];
  c_expl_temp.colidx->size[0] = e_expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(c_expl_temp.colidx, i);
  loop_ub = e_expl_temp.colidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    c_expl_temp.colidx->data[i] = e_expl_temp.colidx->data[i];
  }

  i = c_expl_temp.rowidx->size[0];
  c_expl_temp.rowidx->size[0] = e_expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(c_expl_temp.rowidx, i);
  loop_ub = e_expl_temp.rowidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    c_expl_temp.rowidx->data[i] = e_expl_temp.rowidx->data[i];
  }

  i = b_tempd->size[0] * b_tempd->size[1];
  b_tempd->size[0] = 1;
  b_tempd->size[1] = (int)N;
  emxEnsureCapacity_real_T(b_tempd, i);
  loop_ub = (int)N;
  for (i = 0; i < loop_ub; i++) {
    b_tempd->data[i] = 1.0;
  }

  sparse(K, l * N, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  g_diag(K,q, dv);
  f_kron(K,b_tempd, dv, r2);
  p_sparse_horzcat(K,b_i, t2_n, c_expl_temp.d, c_expl_temp.colidx,
                   c_expl_temp.rowidx, e_expl_temp.m, e_expl_temp.n, r2,
                   input_sizes_idx_1, b_sizes_idx_1, &A102);
  loop_ub = A102.d->size[0];
  emxFree_real_T(&r2);
  for (i = 0; i < loop_ub; i++) {
    A102.d->data[i] = -A102.d->data[i];
  }

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* Aeq6 */
  /* \sum_{k=1}^K\ S_k^n= 0, for n \in t^{pi} */
  i = templ->size[0] * templ->size[1];
  templ->size[0] = (int)tpi_size;
  templ->size[1] = i1;
  emxEnsureCapacity_real_T(templ, i);
  loop_ub = (int)tpi_size * i1;
  for (i = 0; i < loop_ub; i++) {
    templ->data[i] = 0.0;
  }

  /*  for S */
  i = (int)tpi_size;
  for (b_i = 0; b_i < i; b_i++) {
    b_appo = (3.0 * l * N + rt_powd_snf(l, 2.0) * N * M) + l * N * M;
    b_y = tpi[b_i];
    b_a = (b_appo + (b_y - 1.0) * K) + 1.0;
    b_appo += b_y * K;
    if (b_a > b_appo) {
      i1 = 0;
      Ay_n = 0;
    } else {
      i1 = (int)b_a - 1;
      Ay_n = (int)b_appo;
    }

    loop_ub = Ay_n - i1;
    for (Ay_n = 0; Ay_n < loop_ub; Ay_n++) {
      templ->data[b_i + templ->size[0] * (i1 + Ay_n)] = 1.0;
    }
  }

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* Aeq7 */
  /* \sum_{k\in\K_i}\sum_{n=1}^{ceil{\tau_i/\delta t}-1} S_k^n=0 \for i=1,\dots,l */
  i = temp->size[0] * temp->size[1];
  temp->size[0] = (int)l;
  temp->size[1] = i2;
  emxEnsureCapacity_real_T(temp, i);
  loop_ub = (int)l * i2;
  for (i = 0; i < loop_ub; i++) {
    temp->data[i] = 0.0;
  }

  i = (int)l;
  for (b_i = 0; b_i < i; b_i++) {
    i1 = (int)K;
    for (sizes_idx_1 = 0; sizes_idx_1 < i1; sizes_idx_1++) {
      if (Ki[sizes_idx_1] == (double)b_i + 1.0) {
        i2 = (int)(delay[b_i] - 1.0);
        for (input_sizes_idx_1 = 0; input_sizes_idx_1 < i2; input_sizes_idx_1++)
        {
          temp->data[b_i + temp->size[0] * ((int)(((double)sizes_idx_1 + 1.0) +
            (((double)input_sizes_idx_1 + 1.0) - 1.0) * K) - 1)] = 1.0;
        }
      }
    }
  }

  emxInitStruct_sparse(&Aeq7);
  sparse(l, (3.0 * l * N + N * M * l) + N * M * l * l, e_expl_temp.d,
         e_expl_temp.colidx, e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  sparse(l, K * ((N - maxval_tmp) + 1.0), d_expl_temp.d, d_expl_temp.colidx,
         d_expl_temp.rowidx, &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  sparse(l, K * N, c_expl_temp.d, c_expl_temp.colidx, c_expl_temp.rowidx, &Ay_m,
         &Ay_n, &c_i);
  sparse(l, l * N, b_expl_temp.d, b_expl_temp.colidx, b_expl_temp.rowidx, &k,
         &sizes_idx_1, &c_i);
  k_sparse_horzcat(b_i, t2_n, temp, d_expl_temp.d, d_expl_temp.colidx,
                   d_expl_temp.rowidx, input_sizes_idx_1, b_sizes_idx_1, Ay_m,
                   Ay_n, k, sizes_idx_1, &Aeq7);

  /* %%%%%%%% */
  /* Az1  -z_k-\sum_{n=1}^N n S_k^n\le -s_k  \forall k=1,\dots,K\\ */
  if (N < 1.0) {
    b_tempd->size[0] = 1;
    b_tempd->size[1] = 0;
  } else if (rtIsInf(N) && (1.0 == N)) {
    i = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = 1;
    emxEnsureCapacity_real_T(b_tempd, i);
    b_tempd->data[0] = rtNaN;
  } else {
    i = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = (int)floor(N - 1.0) + 1;
    emxEnsureCapacity_real_T(b_tempd, i);
    loop_ub = (int)floor(N - 1.0);
    for (i = 0; i <= loop_ub; i++) {
      b_tempd->data[i] = (double)i + 1.0;
    }
  }

  sparse(K, (3.0 * l * N + N * M * l * l) + l * M * N, e_expl_temp.d,
         e_expl_temp.colidx, e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  eye(K, r1);
  e_kron(b_tempd, r1, b_r);
  sparse(K, K * N + l * N, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  eye(K, temp);
  i = b_delay_mat->size[0] * b_delay_mat->size[1];
  b_delay_mat->size[0] = b_r->size[0];
  b_delay_mat->size[1] = b_r->size[1];
  emxEnsureCapacity_real_T(b_delay_mat, i);
  loop_ub = b_r->size[0] * b_r->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_delay_mat->data[i] = -b_r->data[i];
  }

  i = b_temp->size[0] * b_temp->size[1];
  b_temp->size[0] = temp->size[0];
  b_temp->size[1] = temp->size[1];
  emxEnsureCapacity_real_T(b_temp, i);
  loop_ub = temp->size[0] * temp->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_temp->data[i] = -temp->data[i];
  }

  q_sparse_horzcat(b_i, t2_n, b_delay_mat, input_sizes_idx_1, b_sizes_idx_1,
                   b_temp, &e_expl_temp);
  i = expl_temp.d->size[0];
  expl_temp.d->size[0] = e_expl_temp.d->size[0];
  emxEnsureCapacity_real_T(expl_temp.d, i);
  loop_ub = e_expl_temp.d->size[0];
  emxFree_real_T(&b_delay_mat);
  for (i = 0; i < loop_ub; i++) {
    expl_temp.d->data[i] = e_expl_temp.d->data[i];
  }

  i = expl_temp.colidx->size[0];
  expl_temp.colidx->size[0] = e_expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(expl_temp.colidx, i);
  loop_ub = e_expl_temp.colidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    expl_temp.colidx->data[i] = e_expl_temp.colidx->data[i];
  }

  i = expl_temp.rowidx->size[0];
  expl_temp.rowidx->size[0] = e_expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(expl_temp.rowidx, i);
  loop_ub = e_expl_temp.rowidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    expl_temp.rowidx->data[i] = e_expl_temp.rowidx->data[i];
  }

  aeq1_temp1_m = e_expl_temp.m;
  aeq1_temp1_n = e_expl_temp.n;

  /* Az2  -z_k+\sum_{n=1}^N n S_k^n\le s_k  \forall k=1,\dots,K\\ */
  if (N < 1.0) {
    b_tempd->size[0] = 1;
    b_tempd->size[1] = 0;
  } else if (rtIsInf(N) && (1.0 == N)) {
    i = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = 1;
    emxEnsureCapacity_real_T(b_tempd, i);
    b_tempd->data[0] = rtNaN;
  } else {
    i = b_tempd->size[0] * b_tempd->size[1];
    b_tempd->size[0] = 1;
    b_tempd->size[1] = (int)floor(N - 1.0) + 1;
    emxEnsureCapacity_real_T(b_tempd, i);
    loop_ub = (int)floor(N - 1.0);
    for (i = 0; i <= loop_ub; i++) {
      b_tempd->data[i] = (double)i + 1.0;
    }
  }

  sparse(K, (3.0 * l * N + N * M * l * l) + l * M * N, e_expl_temp.d,
         e_expl_temp.colidx, e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  sparse(K, K * N + l * N, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  eye(K, temp);
  i = b_temp->size[0] * b_temp->size[1];
  b_temp->size[0] = temp->size[0];
  b_temp->size[1] = temp->size[1];
  emxEnsureCapacity_real_T(b_temp, i);
  loop_ub = temp->size[0] * temp->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_temp->data[i] = -temp->data[i];
  }

  emxFree_real_T(&temp);
  eye(K, b_r);
  e_kron(b_tempd, b_r, r1);
  q_sparse_horzcat(b_i, t2_n, r1, input_sizes_idx_1, b_sizes_idx_1, b_temp,
                   &aeq1_temp2);

  /* Ay  y_k+\sum_{n=1}^N D_k^n=d_k \forall k=1,\dots,K\\ */
  sparse(K, ((3.0 * l * N + N * M * l * l) + l * M * N) + K * N, e_expl_temp.d,
         e_expl_temp.colidx, e_expl_temp.rowidx, &b_i, &t2_n, &c_i);
  i = b_tempd->size[0] * b_tempd->size[1];
  b_tempd->size[0] = 1;
  b_tempd->size[1] = (int)N;
  emxEnsureCapacity_real_T(b_tempd, i);
  loop_ub = (int)N;
  emxFree_real_T(&b_temp);
  for (i = 0; i < loop_ub; i++) {
    b_tempd->data[i] = 1.0;
  }

  sparse(K, K + l * N, d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
         &input_sizes_idx_1, &b_sizes_idx_1, &c_i);
  eye(K, b_r);
  e_kron(b_tempd, b_r, r1);
  eye(K, b_r);
  q_sparse_horzcat(b_i, t2_n, r1, input_sizes_idx_1, b_sizes_idx_1, b_r,
                   &e_expl_temp);
  i = c_expl_temp.d->size[0];
  c_expl_temp.d->size[0] = e_expl_temp.d->size[0];
  emxEnsureCapacity_real_T(c_expl_temp.d, i);
  loop_ub = e_expl_temp.d->size[0];
  emxFree_real_T(&r1);
  emxFree_real_T(&b_r);
  emxFree_real_T(&b_tempd);
  for (i = 0; i < loop_ub; i++) {
    c_expl_temp.d->data[i] = e_expl_temp.d->data[i];
  }

  i = c_expl_temp.colidx->size[0];
  c_expl_temp.colidx->size[0] = e_expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(c_expl_temp.colidx, i);
  loop_ub = e_expl_temp.colidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    c_expl_temp.colidx->data[i] = e_expl_temp.colidx->data[i];
  }

  i = c_expl_temp.rowidx->size[0];
  c_expl_temp.rowidx->size[0] = e_expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(c_expl_temp.rowidx, i);
  loop_ub = e_expl_temp.rowidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    c_expl_temp.rowidx->data[i] = e_expl_temp.rowidx->data[i];
  }

  Ay_m = e_expl_temp.m;
  Ay_n = e_expl_temp.n;

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* A8=[sparse(N*l,N*l),eye(N*l),sparse(N*l,2*N*l+l*l*N*M+l*N*M+2*K*N)]; */
  /* b8=repmat(ubL,[N,1]); */
  loop_ub = (int)(((a + K) + K) + varargin_1_tmp);
  i = lb->size[0];
  lb->size[0] = loop_ub;
  emxEnsureCapacity_real_T(lb, i);
  for (i = 0; i < loop_ub; i++) {
    lb->data[i] = 0.0;
  }

  maxval_tmp = b_maximum(N,r);
  loop_ub = (int)(2.0 * l * N);
  i = b->size[0];
  b->size[0] = loop_ub;
  emxEnsureCapacity_real_T(b, i);
  for (i = 0; i < loop_ub; i++) {
    b->data[i] = 1.0;
  }

  loop_ub = (int)((((varargin_1_tmp + y) + b_d) + d1) + varargin_1_tmp);
  i = ub->size[0];
  ub->size[0] = ((b->size[0] + loop_ub) + (int)K) + K;
  emxEnsureCapacity_real_T(ub, i);
  sizes_idx_1 = b->size[0];
  for (i = 0; i < sizes_idx_1; i++) {
    ub->data[i] = maxval_tmp * b->data[i];
  }

  for (i = 0; i < loop_ub; i++) {
    ub->data[i + b->size[0]] = 1.0;
  }

  sizes_idx_1 = (int)K;
  for (i = 0; i < sizes_idx_1; i++) {
    ub->data[(i + b->size[0]) + loop_ub] = N;
  }

  for (i = 0; i < K; i++) {
    ub->data[((i + b->size[0]) + loop_ub) + (int)K] = d[i];
  }

  /*  Preset originale */
  /*  A=[A11;A12;A2;A3;A4;A6;A71;A91;A92;A22;A8]; */
  /*  b=[b11;b12;b2;b3;b4;b6;b71;b91;b92;b22;b8]; */
  /*  Aeq=[Aeq1;Aeq2;Aeq3;Aeq4;Aeq5;Aeq6;Aeq7;A2eqini]; */
  /*  beq=[beq1;beq2;beq3;beq4;beq5;beq6;beq7;b2eqini]; */
  /* %%% */
  c_sparse_vertcat(g_expl_temp.d, g_expl_temp.colidx, g_expl_temp.rowidx, A11_m,
                   A11_n, h_expl_temp.d, h_expl_temp.colidx, h_expl_temp.rowidx,
                   A12_m, A12_n, i_expl_temp.d, i_expl_temp.colidx,
                   i_expl_temp.rowidx, i_expl_temp.m, i_expl_temp.n, A22, A3, A4,
                   A5, A6, A7, Aeq1, Aeqini, Aeqini1, A102, &b_expl_temp);
  d_sparse_vertcat(f_expl_temp.d, f_expl_temp.colidx, f_expl_temp.rowidx,
                   f_expl_temp.m, f_expl_temp.n, Aeq2, j_expl_temp.d,
                   j_expl_temp.colidx, j_expl_temp.rowidx, j_expl_temp.m,
                   j_expl_temp.n, k_expl_temp.d, k_expl_temp.colidx,
                   k_expl_temp.rowidx, k_expl_temp.m, k_expl_temp.n, Aeq5, templ,
                   Aeq7, A2eqini, &Aeq1);

  /* %%% */
  /* NN=4*l*N+l^2*N*M+l*N*M+2*K*N; */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  sparse(b_expl_temp.m, K, e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx,
         &b_i, &t2_n, &c_i);
  r_sparse_horzcat(b_expl_temp.d, b_expl_temp.colidx, b_expl_temp.rowidx,
                   b_expl_temp.m, b_expl_temp.n, b_i, t2_n, &d_expl_temp);
  sparse_vertcat(expl_temp.d, expl_temp.colidx, expl_temp.rowidx, aeq1_temp1_m,
                 aeq1_temp1_n, aeq1_temp2.d, aeq1_temp2.colidx,
                 aeq1_temp2.rowidx, aeq1_temp2.m, aeq1_temp2.n, &e_expl_temp);
  sparse_vertcat(d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
                 d_expl_temp.m, d_expl_temp.n, e_expl_temp.d, e_expl_temp.colidx,
                 e_expl_temp.rowidx, e_expl_temp.m, e_expl_temp.n, &f_expl_temp);
  i = b_expl_temp.d->size[0];
  b_expl_temp.d->size[0] = f_expl_temp.d->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.d, i);
  loop_ub = f_expl_temp.d->size[0];
  emxFreeStruct_sparse(&k_expl_temp);
  emxFreeStruct_sparse(&j_expl_temp);
  emxFreeStruct_sparse(&i_expl_temp);
  emxFreeStruct_sparse(&h_expl_temp);
  emxFreeStruct_sparse(&g_expl_temp);
  emxFreeStruct_sparse(&Aeq7);
  emxFreeStruct_sparse(&A102);
  emxFreeStruct_sparse(&Aeq5);
  emxFreeStruct_sparse(&A7);
  emxFreeStruct_sparse(&A6);
  emxFreeStruct_sparse(&A5);
  emxFreeStruct_sparse(&A4);
  emxFreeStruct_sparse(&A3);
  emxFree_real_T(&Aeq2);
  emxFreeStruct_sparse(&A2eqini);
  emxFreeStruct_sparse(&A22);
  emxFreeStruct_sparse(&Aeqini1);
  emxFreeStruct_sparse(&Aeqini);
  emxFree_real_T(&templ);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.d->data[i] = f_expl_temp.d->data[i];
  }

  i = b_expl_temp.colidx->size[0];
  b_expl_temp.colidx->size[0] = f_expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(b_expl_temp.colidx, i);
  loop_ub = f_expl_temp.colidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.colidx->data[i] = f_expl_temp.colidx->data[i];
  }

  i = b_expl_temp.rowidx->size[0];
  b_expl_temp.rowidx->size[0] = f_expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(b_expl_temp.rowidx, i);
  loop_ub = f_expl_temp.rowidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.rowidx->data[i] = f_expl_temp.rowidx->data[i];
  }

  sparse(Aeq1.m, K, e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx, &b_i,
         &t2_n, &c_i);
  r_sparse_horzcat(Aeq1.d, Aeq1.colidx, Aeq1.rowidx, Aeq1.m, Aeq1.n, b_i, t2_n,
                   &e_expl_temp);
  i = expl_temp.d->size[0];
  expl_temp.d->size[0] = e_expl_temp.d->size[0];
  emxEnsureCapacity_real_T(expl_temp.d, i);
  loop_ub = e_expl_temp.d->size[0];
  for (i = 0; i < loop_ub; i++) {
    expl_temp.d->data[i] = e_expl_temp.d->data[i];
  }

  i = expl_temp.colidx->size[0];
  expl_temp.colidx->size[0] = e_expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(expl_temp.colidx, i);
  loop_ub = e_expl_temp.colidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    expl_temp.colidx->data[i] = e_expl_temp.colidx->data[i];
  }

  i = expl_temp.rowidx->size[0];
  expl_temp.rowidx->size[0] = e_expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(expl_temp.rowidx, i);
  loop_ub = e_expl_temp.rowidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    expl_temp.rowidx->data[i] = e_expl_temp.rowidx->data[i];
  }

  /* NN=NN+K; */
  sparse(f_expl_temp.m, K, e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx,
         &b_i, &t2_n, &c_i);
  r_sparse_horzcat(b_expl_temp.d, b_expl_temp.colidx, b_expl_temp.rowidx,
                   f_expl_temp.m, f_expl_temp.n, b_i, t2_n, &aeq1_temp2);
  sparse(e_expl_temp.m, K, e_expl_temp.d, e_expl_temp.colidx, e_expl_temp.rowidx,
         &b_i, &t2_n, &c_i);
  r_sparse_horzcat(expl_temp.d, expl_temp.colidx, expl_temp.rowidx,
                   e_expl_temp.m, e_expl_temp.n, b_i, t2_n, &e_expl_temp);
  i = d_expl_temp.d->size[0];
  d_expl_temp.d->size[0] = e_expl_temp.d->size[0];
  emxEnsureCapacity_real_T(d_expl_temp.d, i);
  loop_ub = e_expl_temp.d->size[0];
  emxFreeStruct_sparse(&expl_temp);
  emxFreeStruct_sparse(&f_expl_temp);
  for (i = 0; i < loop_ub; i++) {
    d_expl_temp.d->data[i] = e_expl_temp.d->data[i];
  }

  i = d_expl_temp.colidx->size[0];
  d_expl_temp.colidx->size[0] = e_expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(d_expl_temp.colidx, i);
  loop_ub = e_expl_temp.colidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    d_expl_temp.colidx->data[i] = e_expl_temp.colidx->data[i];
  }

  i = d_expl_temp.rowidx->size[0];
  d_expl_temp.rowidx->size[0] = e_expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(d_expl_temp.rowidx, i);
  loop_ub = e_expl_temp.rowidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    d_expl_temp.rowidx->data[i] = e_expl_temp.rowidx->data[i];
  }

  sparse_vertcat(d_expl_temp.d, d_expl_temp.colidx, d_expl_temp.rowidx,
                 e_expl_temp.m, e_expl_temp.n, c_expl_temp.d, c_expl_temp.colidx,
                 c_expl_temp.rowidx, Ay_m, Ay_n, &Aeq1);

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  b_appo = ((j_funct[0] + j_funct[1]) + j_funct[2]) + j_funct[3];

  /* f=[zeros(l*N,1); j(3)*1/(sum(r)*l*N)*ones(l*N,1); zeros(l*N,1);j(4)*kron(ones(N*M,1),psi_temp)/(PSI*l*l*M*N); zeros(l*N*M,1) ; zeros(N*K,1) ; zeros(N*K,1);zeros(N*l,1);j(1)*alfa.*ones(K,1)/(sum(Dt)*K);j(2)*beta.*ones(K,1)/(sum(Dv)*K)]; */
  y = r[0];
  emxFreeStruct_sparse(&d_expl_temp);
  emxFreeStruct_sparse(&c_expl_temp);
  for (k = 0; k < N-1; k++) {
    y += r[k + 1];
  }

  maxval_tmp = j_funct[2] / b_appo / y;
  i = b->size[0];
  b->size[0] = (int)varargin_1_tmp;
  emxEnsureCapacity_real_T(b, i);
  for (i = 0; i < index_canc_idx_0; i++) {
    b->data[i] = 1.0;
  }

  a = j_funct[3] / b_appo;
  k = (int)d2;
  i = b_templ->size[0];
  b_templ->size[0] = (int)d2;
  emxEnsureCapacity_real_T(b_templ, i);
  for (i = 0; i < k; i++) {
    b_templ->data[i] = 1.0;
  }

  g_kron(l,b_templ, psi, t4_d);
  b_a = j_funct[0] / b_appo;
  b_appo = j_funct[1] / b_appo;
  y = Dt[0];
  b_y = Dv[0];
  for (k = 0; k < 18; k++) {
    y += Dt[k + 1];
    b_y += Dv[k + 1];
  }

  i = f->size[0];
  loop_ub = (int)(N * K);
  f->size[0] = ((((((((int)varargin_1_tmp + b->size[0]) + (int)varargin_1_tmp) +
                    t4_d->size[0]) + loop_ub_tmp) + (int)varargin_1) + loop_ub)
                + (int)varargin_1_tmp) + 2*K;
  emxEnsureCapacity_real_T(f, i);
  for (i = 0; i < index_canc_idx_0; i++) {
    f->data[i] = 0.0;
  }

  sizes_idx_1 = b->size[0];
  for (i = 0; i < sizes_idx_1; i++) {
    f->data[i + (int)varargin_1_tmp] = maxval_tmp * b->data[i];
  }

  for (i = 0; i < index_canc_idx_0; i++) {
    f->data[(i + (int)varargin_1_tmp) + b->size[0]] = 0.0;
  }

  sizes_idx_1 = t4_d->size[0];
  for (i = 0; i < sizes_idx_1; i++) {
    f->data[((i + (int)varargin_1_tmp) + b->size[0]) + (int)varargin_1_tmp] = a *
      t4_d->data[i] / N;
  }

  for (i = 0; i < loop_ub_tmp; i++) {
    f->data[(((i + (int)varargin_1_tmp) + b->size[0]) + (int)varargin_1_tmp) +
      t4_d->size[0]] = 0.0;
  }

  for (i = 0; i < loop_ub; i++) {
    f->data[((((i + (int)varargin_1_tmp) + b->size[0]) + (int)varargin_1_tmp) +
             t4_d->size[0]) + loop_ub_tmp] = 0.0;
  }

  for (i = 0; i < loop_ub; i++) {
    f->data[(((((i + (int)varargin_1_tmp) + b->size[0]) + (int)varargin_1_tmp) +
              t4_d->size[0]) + loop_ub_tmp) + loop_ub] = 0.0;
  }

  for (i = 0; i < index_canc_idx_0; i++) {
    f->data[((((((i + (int)varargin_1_tmp) + b->size[0]) + (int)varargin_1_tmp)
               + t4_d->size[0]) + loop_ub_tmp) + loop_ub) + loop_ub] = 0.0;
  }

  for (i = 0; i < K; i++) {
    f->data[(((((((i + (int)varargin_1_tmp) + b->size[0]) + (int)varargin_1_tmp)
                + t4_d->size[0]) + loop_ub_tmp) + loop_ub) + loop_ub) + (int)
      varargin_1_tmp] = b_a * alfa[i] / y;
  }

  for (i = 0; i < K; i++) {
    f->data[((((((((i + (int)varargin_1_tmp) + b->size[0]) + (int)varargin_1_tmp)
                 + t4_d->size[0]) + loop_ub_tmp) + loop_ub) + loop_ub) + (int)
             varargin_1_tmp) + K] = b_appo * beta[i] / b_y;
  }

  emxInit_char_T(&ctypePL, 2);

  /* f=[zeros(l*N,1); j2*  1/sum(r)  *ones(l*N,1); zeros(l*N,1);j3*  kron(ones(N*M,1),psi_temp)/PSI  ; zeros(l*N*M,1) ; zeros(N*K,1) ; zeros(N*K,1);zeros(N*l,1);j1  *alfa.*ones(K,1)/ sum(Dt) ;j11 *beta.*ones(K,1)/sum(Dv)]; */
  i = ctypePL->size[0] * ctypePL->size[1];
  ctypePL->size[0] = 1;
  ctypePL->size[1] = (int)varargin_1_tmp;
  emxEnsureCapacity_char_T(ctypePL, i);
  k = -1;
  i = (int)varargin_1_tmp;
  for (input_sizes_idx_1 = 0; input_sizes_idx_1 < i; input_sizes_idx_1++) {
    k++;
    ctypePL->data[k] = 'C';
  }

  emxInit_char_T(&ctypeE, 2);
  i = ctypeE->size[0] * ctypeE->size[1];
  ctypeE->size[0] = 1;
  ctypeE->size[1] = (int)index_canc_idx_0_tmp;
  emxEnsureCapacity_char_T(ctypeE, i);
  k = -1;
  i = (int)index_canc_idx_0_tmp;
  for (input_sizes_idx_1 = 0; input_sizes_idx_1 < i; input_sizes_idx_1++) {
    k++;
    ctypeE->data[k] = 'B';
  }

  emxInit_char_T(&ctypeS, 2);
  i = ctypeS->size[0] * ctypeS->size[1];
  ctypeS->size[0] = 1;
  ctypeS->size[1] = (int)varargin_1;
  emxEnsureCapacity_char_T(ctypeS, i);
  k = -1;
  i = (int)varargin_1;
  for (input_sizes_idx_1 = 0; input_sizes_idx_1 < i; input_sizes_idx_1++) {
    k++;
    ctypeS->data[k] = 'B';
  }

  emxInit_char_T(&ctypeZ, 2);
  i = ctypeZ->size[0] * ctypeZ->size[1];
  ctypeZ->size[0] = 1;
  ctypeZ->size[1] = (int)K;
  emxEnsureCapacity_char_T(ctypeZ, i);
  k = -1;
  i = (int)K;
  for (input_sizes_idx_1 = 0; input_sizes_idx_1 < i; input_sizes_idx_1++) {
    k++;
    ctypeZ->data[k] = 'C';
  }

  i = (int)tp_size;
  if (0 <= (int)tp_size - 1) {
    c_varargin_1 = b_varargin_1;
    i3 = (int)b_varargin_1;
  }

  emxInit_char_T(&c_y, 2);
  for (b_i = 0; b_i < i; b_i++) {
    b_d = tp[b_i];
    d1 = l * (b_d - 1.0) * M + 1.0;
    if (d1 > l * b_d * M) {
      i1 = 1;
    } else {
      i1 = (int)d1;
    }

    i2 = c_y->size[0] * c_y->size[1];
    c_y->size[0] = 1;
    c_y->size[1] = (int)c_varargin_1;
    emxEnsureCapacity_char_T(c_y, i2);
    k = -1;
    for (input_sizes_idx_1 = 0; input_sizes_idx_1 < i3; input_sizes_idx_1++) {
      k++;
      c_y->data[k] = 'C';
    }

    loop_ub = c_y->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      ctypeE->data[(i1 + i2) - 1] = c_y->data[i2];
    }
  }

  /* ctype=strcat(ctypePL,ctypePL,ctypeG,ctypeF,ctypeE,ctypeS,ctypeD,ctypePL,ctypeZ,ctypeZ); */
  varargin_1 = l * l * M * N;
  i = c_y->size[0] * c_y->size[1];
  c_y->size[0] = 1;
  c_y->size[1] = (int)varargin_1;
  emxEnsureCapacity_char_T(c_y, i);
  k = -1;
  i = (int)varargin_1;
  for (input_sizes_idx_1 = 0; input_sizes_idx_1 < i; input_sizes_idx_1++) {
    k++;
    c_y->data[k] = 'C';
  }

  b_sprintf(ctypePL, ctypePL, ctypePL, c_y, ctypeE, ctypeS, ctypeS, ctypePL,
            ctypeZ, ctypeZ, ctype);
  c_i = aeq1_temp2.m;
  repelem(aeq1_temp2.m, ctypePL);
  b_repelem(Aeq1.m, ctypeE);
  c_sprintf(ctypePL, ctypeE, A_sense);
  sparse_vertcat(aeq1_temp2.d, aeq1_temp2.colidx, aeq1_temp2.rowidx,
                 aeq1_temp2.m, aeq1_temp2.n, Aeq1.d, Aeq1.colidx, Aeq1.rowidx,
                 Aeq1.m, Aeq1.n, &e_expl_temp);
  i = b_expl_temp.d->size[0];
  b_expl_temp.d->size[0] = e_expl_temp.d->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.d, i);
  loop_ub = e_expl_temp.d->size[0];
  emxFree_char_T(&c_y);
  emxFree_char_T(&ctypeZ);
  emxFree_char_T(&ctypeS);
  emxFree_char_T(&ctypeE);
  emxFree_char_T(&ctypePL);
  emxFreeStruct_sparse(&Aeq1);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.d->data[i] = e_expl_temp.d->data[i];
  }

  i = b_expl_temp.colidx->size[0];
  b_expl_temp.colidx->size[0] = e_expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(b_expl_temp.colidx, i);
  loop_ub = e_expl_temp.colidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.colidx->data[i] = e_expl_temp.colidx->data[i];
  }

  i = b_expl_temp.rowidx->size[0];
  b_expl_temp.rowidx->size[0] = e_expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(b_expl_temp.rowidx, i);
  loop_ub = e_expl_temp.rowidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.rowidx->data[i] = e_expl_temp.rowidx->data[i];
  }

  i = b->size[0];
  b->size[0] = (((((((((((((((((((appo_rep->size[0] + beqini->size[0]) +
    b21->size[0]) + b22->size[0]) + (int)M) + b4->size[0]) + (int)M) + (int)N) +
    (int)K) + b8->size[0]) + (int)K) + (int)K) + beq1->size[0]) + (int)tp_size)
                     + beq3->size[0]) + beq4->size[0]) + (int)K) + (int)tpi_size)
                 + (int)l) + (int)N) + 2*l+4*K;
  emxEnsureCapacity_real_T(b, i);
  for (i = 0; i < l; i++) {
    b->data[i] = H0[i];
  }

  loop_ub = appo_rep->size[0];
  for (i = 0; i < loop_ub; i++) {
    b->data[i + l] = appo_rep->data[i] * eps_min;
  }

  for (i = 0; i < l; i++) {
    b->data[(i + appo_rep->size[0]) + l] = -H0[i];
  }

  loop_ub = beqini->size[0];
  for (i = 0; i < loop_ub; i++) {
    b->data[(i + appo_rep->size[0]) + 2*l] = beqini->data[i] * eps_min;
  }

  loop_ub = b21->size[0];
  for (i = 0; i < loop_ub; i++) {
    b->data[((i + appo_rep->size[0]) + beqini->size[0]) + 2*l] = b21->data[i];
  }

  loop_ub = b22->size[0];
  for (i = 0; i < loop_ub; i++) {
    b->data[(((i + appo_rep->size[0]) + beqini->size[0]) + b21->size[0]) + 2*l] =
      b22->data[i];
  }

  loop_ub = (int)M;
  for (i = 0; i < loop_ub; i++) {
    b->data[((((i + appo_rep->size[0]) + beqini->size[0]) + b21->size[0]) +
             b22->size[0]) + 2*l] = 1.0;
  }

  loop_ub = b4->size[0];
  for (i = 0; i < loop_ub; i++) {
    b->data[(((((i + appo_rep->size[0]) + beqini->size[0]) + b21->size[0]) +
              b22->size[0]) + (int)M) + 2*l] = b4->data[i];
  }

  loop_ub = (int)M;
  for (i = 0; i < loop_ub; i++) {
    b->data[((((((i + appo_rep->size[0]) + beqini->size[0]) + b21->size[0]) +
               b22->size[0]) + (int)M) + b4->size[0]) + 2*l] = 0.0;
  }

  loop_ub = (int)N;
  for (i = 0; i < loop_ub; i++) {
    b->data[(((((((i + appo_rep->size[0]) + beqini->size[0]) + b21->size[0]) +
                b22->size[0]) + (int)M) + b4->size[0]) + (int)M) + 2*l] = 1.0;
  }

  loop_ub = (int)K;
  for (i = 0; i < loop_ub; i++) {
    b->data[((((((((i + appo_rep->size[0]) + beqini->size[0]) + b21->size[0]) +
                 b22->size[0]) + (int)M) + b4->size[0]) + (int)M) + (int)N) + 2*l]
      = 1.0;
  }

  loop_ub = b8->size[0];
  for (i = 0; i < loop_ub; i++) {
    b->data[(((((((((i + appo_rep->size[0]) + beqini->size[0]) + b21->size[0]) +
                  b22->size[0]) + (int)M) + b4->size[0]) + (int)M) + (int)N) +
             (int)K) + 2*l] = b8->data[i];
  }

  loop_ub = (int)K;
  for (i = 0; i < loop_ub; i++) {
    b->data[((((((((((i + appo_rep->size[0]) + beqini->size[0]) + b21->size[0])
                   + b22->size[0]) + (int)M) + b4->size[0]) + (int)M) + (int)N)
              + (int)K) + b8->size[0]) + 2*l] = N;
  }

  for (i = 0; i < K; i++) {
    b->data[(((((((((((i + appo_rep->size[0]) + beqini->size[0]) + b21->size[0])
                    + b22->size[0]) + (int)M) + b4->size[0]) + (int)M) + (int)N)
               + (int)K) + b8->size[0]) + (int)K) + 2*l] = q[i] * d[i];
  }

  loop_ub = (int)K;
  for (i = 0; i < loop_ub; i++) {
    b->data[(((((((((((i + appo_rep->size[0]) + beqini->size[0]) + b21->size[0])
                    + b22->size[0]) + (int)M) + b4->size[0]) + (int)M) + (int)N)
               + (int)K) + b8->size[0]) + (int)K) + 2*l+K] = 0.0;
  }

  for (i = 0; i < K; i++) {
    b_d = s[i];
    b->data[(((((((((((i + (appo_rep->size[0] + beqini->size[0])) + b21->size[0])
                     + b22->size[0]) + (int)M) + b4->size[0]) + (int)M) + (int)N)
                + (int)K) + b8->size[0]) + (int)K) + (int)K) + 33] = -b_d;
    b->data[(((((((((((i + (appo_rep->size[0] + beqini->size[0])) + b21->size[0])
                     + b22->size[0]) + (int)M) + b4->size[0]) + (int)M) + (int)N)
                + (int)K) + b8->size[0]) + (int)K) + (int)K) + 2*l+2*K] = b_d;
  }

  loop_ub = beq1->size[0];
  for (i = 0; i < loop_ub; i++) {
    b->data[(((((((((((i + (appo_rep->size[0] + beqini->size[0])) + b21->size[0])
                     + b22->size[0]) + (int)M) + b4->size[0]) + (int)M) + (int)N)
                + (int)K) + b8->size[0]) + (int)K) + (int)K) + 2*l+3*K] = beq1->
      data[i];
  }

  loop_ub = (int)tp_size;
  for (i = 0; i < loop_ub; i++) {
    b->data[(((((((((((i + ((appo_rep->size[0] + beqini->size[0]) + b21->size[0]))
                      + b22->size[0]) + (int)M) + b4->size[0]) + (int)M) + (int)
                  N) + (int)K) + b8->size[0]) + (int)K) + (int)K) + beq1->size[0])
      + 2*l+3*K] = 0.0;
  }

  loop_ub = beq3->size[0];
  for (i = 0; i < loop_ub; i++) {
    b->data[(((((((((((i + (((appo_rep->size[0] + beqini->size[0]) + b21->size[0])
      + b22->size[0])) + (int)M) + b4->size[0]) + (int)M) + (int)N) + (int)K) +
                 b8->size[0]) + (int)K) + (int)K) + beq1->size[0]) + (int)
             tp_size) + 2*l+3*K] = beq3->data[i];
  }

  loop_ub = beq4->size[0];
  for (i = 0; i < loop_ub; i++) {
    b->data[(((((((((((i + ((((appo_rep->size[0] + beqini->size[0]) + b21->size
      [0]) + b22->size[0]) + (int)M)) + b4->size[0]) + (int)M) + (int)N) + (int)
                   K) + b8->size[0]) + (int)K) + (int)K) + beq1->size[0]) + (int)
              tp_size) + beq3->size[0]) + 2*l+3*K] = beq4->data[i];
  }

  loop_ub = (int)K;
  for (i = 0; i < loop_ub; i++) {
    b->data[(((((((((((i + (((((appo_rep->size[0] + beqini->size[0]) + b21->
      size[0]) + b22->size[0]) + (int)M) + b4->size[0])) + (int)M) + (int)N) +
                    (int)K) + b8->size[0]) + (int)K) + (int)K) + beq1->size[0])
               + (int)tp_size) + beq3->size[0]) + beq4->size[0]) + 2*l+3*K] = 0.0;
  }

  loop_ub = (int)tpi_size;
  for (i = 0; i < loop_ub; i++) {
    b->data[(((((((((((i + ((((((appo_rep->size[0] + beqini->size[0]) +
      b21->size[0]) + b22->size[0]) + (int)M) + b4->size[0]) + (int)M)) + (int)N)
                     + (int)K) + b8->size[0]) + (int)K) + (int)K) + beq1->size[0])
                + (int)tp_size) + beq3->size[0]) + beq4->size[0]) + (int)K) + 2*l+3*K]
      = 0.0;
  }

  loop_ub = (int)l;
  for (i = 0; i < loop_ub; i++) {
    b->data[(((((((((((i + (((((((appo_rep->size[0] + beqini->size[0]) +
      b21->size[0]) + b22->size[0]) + (int)M) + b4->size[0]) + (int)M) + (int)N))
                      + (int)K) + b8->size[0]) + (int)K) + (int)K) + beq1->size
                  [0]) + (int)tp_size) + beq3->size[0]) + beq4->size[0]) + (int)
              K) + (int)tpi_size) + 2*l+3*K] = 0.0;
  }

  loop_ub = (int)N;
  for (i = 0; i < loop_ub; i++) {
    b->data[(((((((((((i + ((((((((appo_rep->size[0] + beqini->size[0]) +
      b21->size[0]) + b22->size[0]) + (int)M) + b4->size[0]) + (int)M) + (int)N)
      + (int)K)) + b8->size[0]) + (int)K) + (int)K) + beq1->size[0]) + (int)
                  tp_size) + beq3->size[0]) + beq4->size[0]) + (int)K) + (int)
              tpi_size) + (int)l) + 2*l+3*K] = 0.0;
  }

  for (i = 0; i < K; i++) {
    b->data[(((((((((((i + (((((((((appo_rep->size[0] + beqini->size[0]) +
      b21->size[0]) + b22->size[0]) + (int)M) + b4->size[0]) + (int)M) + (int)N)
      + (int)K) + b8->size[0])) + (int)K) + (int)K) + beq1->size[0]) + (int)
                   tp_size) + beq3->size[0]) + beq4->size[0]) + (int)K) + (int)
               tpi_size) + (int)l) + (int)N) + 2*l+3*K] = d[i];
  }

  emxFree_int8_T(&b8);
  emxFree_int8_T(&beq4);
  emxFree_int8_T(&b4);
  emxFree_int8_T(&beq3);
  emxFree_int8_T(&b22);
  emxFree_int8_T(&b21);
  emxFree_real_T(&beq1);
  emxFree_real_T(&beqini);
  emxFree_real_T(&appo_rep);
  emxInit_boolean_T(&t5_d, 1);
  emxInit_int32_T(&t5_colidx, 1);
  sparse_ne(b_expl_temp.d, b_expl_temp.colidx, b_expl_temp.rowidx, e_expl_temp.m,
            e_expl_temp.n, t5_d, t5_colidx, ini_P_tmp, &sizes_idx_1, &k);
  b_sum(t5_d, t5_colidx, sizes_idx_1, k, t4_d, t4_colidx, t4_rowidx,
        &input_sizes_idx_1);
  sparse_ctranspose(t4_d, t4_colidx, t4_rowidx, input_sizes_idx_1, b_templ,
                    ini_P_tmp, t5_colidx, &sizes_idx_1);
  sparse_full(b_templ, ini_P_tmp, t5_colidx, sizes_idx_1, A_nnz);
  b_eml_find(b_expl_temp.colidx, b_expl_temp.rowidx, ini_P_tmp, t5_colidx);
  i = A_j->size[0];
  A_j->size[0] = t5_colidx->size[0];
  emxEnsureCapacity_real_T(A_j, i);
  loop_ub = t5_colidx->size[0];
  emxFree_real_T(&b_templ);
  emxFreeStruct_sparse(&e_expl_temp);
  emxFree_boolean_T(&t5_d);
  emxFree_int32_T(&t4_rowidx);
  emxFree_int32_T(&t4_colidx);
  emxFree_real_T(&t4_d);
  for (i = 0; i < loop_ub; i++) {
    A_j->data[i] = t5_colidx->data[i];
  }

  emxFree_int32_T(&t5_colidx);
  i = A_indx->size[0];
  A_indx->size[0] = ini_P_tmp->size[0];
  emxEnsureCapacity_real_T(A_indx, i);
  loop_ub = ini_P_tmp->size[0];
  for (i = 0; i < loop_ub; i++) {
    A_indx->data[i] = (double)ini_P_tmp->data[i] - 1.0;
  }

  emxFree_int32_T(&ini_P_tmp);
  sparse_nonzeros(b_expl_temp.d, b_expl_temp.colidx, A_val);
  i = A_beg->size[0];
  A_beg->size[0] = aeq1_temp2.n;
  emxEnsureCapacity_real_T(A_beg, i);
  b_appo = 0.0;
  i = aeq1_temp2.n;
  emxFreeStruct_sparse(&b_expl_temp);
  emxFreeStruct_sparse(&aeq1_temp2);
  for (b_i = 0; b_i < i; b_i++) {
    A_beg->data[b_i] = b_appo;
    b_appo += A_nnz->data[b_i];
  }

  *A_r = c_i;
}

/* End of code generation (get_MILP.c) */
