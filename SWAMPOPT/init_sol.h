//
//  init_sol.h
//  initial_solution
//
//  Created by Valyria on 25/09/20.
//  Copyright © 2020 Valyria. All rights reserved.
//

#ifndef init_sol_h
#define init_sol_h
/* Include Files */
#include <stddef.h>
#include <stdlib.h>

/* Function Declarations */
double* initial_solution(const int N, const int M, const int K, const int l,const int tp_size,const int tpi_size, const int                                 tp[tp_size], const int tpi[tpi_size], const double c[l], const double b_gamma[l], const double rho[l], const double psi[l][l],const int delay[l], const int Ki[K],const int Ii[l], const double q[K], const int s[K], const int d[K], const double eps[K], const double eps_min, int  total_delay,const double H0[l], const double V0[total_delay], const double R0[l]);

#endif /* init_sol_h */
