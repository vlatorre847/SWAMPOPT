//
//  my_functions.h
//  initial_solution
//
//  Created by Valyria on 25/09/20.
//  Copyright © 2020 Valyria. All rights reserved.
//

#ifndef my_functions_h
#define my_functions_h

#include <stdio.h>
#include <stdbool.h>
void check_mat(int n_rows, int n_cols, double mat[n_rows][n_cols]);
void check_p(int n_rows, int n_cols, double **mat);
bool my_ismember(int n,double a,int arr[n]);
double** chan_req(int N,int l,int K,int tp_size,double gamma[l],double q[K],int prec[l],int Ki[K],int tp[tp_size],int delay[l],int** D);
void add_irrigation(int N,int M,int l,int K,int tp_size,double** L,double** P,int** G,double** H,int** D,int stopj[l],int startj[l],double q[K],double c[l],int Ii[l],int Ki[K],double gamma[l],double psi[l][l],int tp[tp_size],double rho[l],int delay[l],const double eps_min,int total_delay,double V0[total_delay],double R0[l]);
int* my_find(int n,int val, int arr[n], int ops, int* size_arg);
int* find_greater(int n,int val, int arr[n], int ops, int* size_arg);
void check_int(int n_rows, int n_cols, int **mat);
int* read_int_vec(char* token, int v_length);
double* read_double_vec(char* token, int v_length);
#endif /* my_functions_h */
