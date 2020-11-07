/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * get_MILP_initialize.c
 *
 * Code generation for function 'get_MILP_initialize'
 *
 */

/* Include files */
#include "get_MILP_initialize.h"
#include "get_MILP_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void get_MILP_initialize(void)
{
  rt_InitInfAndNaN();
  isInitialized_get_MILP = true;
}

/* End of code generation (get_MILP_initialize.c) */
