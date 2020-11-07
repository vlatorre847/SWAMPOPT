//
//  min_max.c
//  initial_solution
//
//  Created by Valyria on 25/09/20.
//  Copyright © 2020 Valyria. All rights reserved.
//

#include "my_math.h"
#include <stdio.h>
#include <stdlib.h>

double my_max(double a,double b){
    if (a>=b)
        return a;
        else
        return b;
}

double my_min(double a,double b){
    if (a<=b)
        return a;
    else
        return b;
}


double my_sum(int n, double a[n]){
    double sum=0;
    for (int i=0;i<n;i++){
        sum=sum+a[i];
        //printf("%f\n",a[i]);
    }
    return sum;
}

int sum_i(int n, int a[n]){
    int sum=0;
    for (int i=0;i<n;i++){
        sum=sum+a[i];
        //printf("%f\n",a[i]);
    }
    return sum;
}
void swap(int* xp, int* yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

// Function to perform Selection Sort
int* my_sort(int arr[], int n, int options)
{
    int i, j, min_idx;
    int* idx_arr=(int*)malloc(n * sizeof(int*));
    for (int i=0;i<n;i++)
    {
        idx_arr[i]=i;
    }
    // One by one move boundary of unsorted subarray
    for (i = 0; i < n - 1; i++) {
        
        // Find the minimum element in unsorted array
        min_idx = i;
        for (j = i + 1; j < n; j++)
            if (options==1)
            {
                if (arr[j] > arr[min_idx])
                    min_idx = j;
            }
            else{
                if (arr[j] < arr[min_idx])
                    min_idx = j;
            }
            if (arr[j] < arr[min_idx])
                min_idx = j;
        
        // Swap the found minimum element
        // with the first element
        swap(&arr[min_idx], &arr[i]);
        swap(&idx_arr[min_idx], &idx_arr[i]);
    }
    return idx_arr;
}
int my_nnz(int n_rows, int n_cols, double** arr){
    int res=0;
    for (int i=0;i<n_rows;i++){
        for (int j=0;j<n_cols;j++){
            if (arr[i][j]!=0)
                res=res+1;
        }
    }
    return res;
}
int nnz_int(int n_rows, int n_cols, int** arr){
    int res=0;
    for (int i=0;i<n_rows;i++){
        for (int j=0;j<n_cols;j++){
            if (arr[i][j]!=0)
                res=res+1;
        }
    }
    return res;
}
