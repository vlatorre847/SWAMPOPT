//
//  min_max.h
//  initial_solution
//
//  Created by Valyria on 25/09/20.
//  Copyright © 2020 Valyria. All rights reserved.
//

#ifndef my_math_h
#define my_math_h

#include <stdio.h>

double my_max(double a,double b);
double my_sum(int n, double a[n]);
double my_min(double a,double b);
int sum_i(int n, int a[n]);
int* my_sort(int arr[], int n, int options);
int my_nnz(int n_rows, int n_cols, double** arr);
int nnz_int(int n_rows, int n_cols, int** arr);
#endif /* min_max_h */
