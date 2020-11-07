
//
//  main.c
//  Irrigazione
//
//  Created by Valyria on 02/09/20.
//  Copyright © 2020 Valyria. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "my_functions.h"

#include "get_MILP.h"
#include "get_MILP_emxAPI.h"
#include "get_MILP_terminate.h"
#include "rt_nonfinite.h"
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "get_MILP_types.h"
#include <ilcplex/cplex.h>
#include "init_sol.h"

int main(int argc, const char * argv[]) {
    FILE* fp1;
    int N;
    int M;
    int l;
    fp1 = fopen("CBEC_2d_conf.txt", "r");
    char temp_l[254];
    fgets(temp_l,254,fp1);
    sscanf(temp_l,"N=%d",&N);
    fgets(temp_l,254,fp1);
    sscanf(temp_l,"M=%d",&M);
    fgets(temp_l,254,fp1);
    sscanf(temp_l,"l=%d",&l);
    
    fgets(temp_l,254,fp1);
    char* token;
    token=strtok(temp_l,"Ii=");
    int* Ii=read_int_vec(token,l);
    
    double dt;
    fgets(temp_l,254,fp1);
    sscanf(temp_l,"dt=%lf",&dt);
    
    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"tau=");
    double* tau=read_double_vec(token,l);;

    double delay[l];
    double temp=0;
    int total_delay=0;
    for (int i = 0; i < l; i++) {
        temp=ceil(tau[i]/dt);
        if (temp<1){
            delay[i]=1;
        }
        else{
            delay[i]=temp;
            total_delay=total_delay+(int)delay[i];
        }
    }
    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"L_in=");
    double* L0=read_double_vec(token,l);

    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"P_in=");
    double* P0=(double*)malloc(total_delay*sizeof(double*));
    int sum_delay=0;
    char* token1=strtok(token,",");
    for(int i=0;i<l;i++){
        double* P_temp=read_double_vec(token1,delay[i]);
        for (int j=0;j<delay[i];j++){
            P0[sum_delay+j]=P_temp[j];
        }
        sum_delay=sum_delay+delay[i];
        token1=strtok(NULL,",");
    }
    for (int i = 0; i < total_delay; i++) {
        P0[i] = P0[i]*60*dt/1000;
    }
    
    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"H_in=");
    double* H0=read_double_vec(token,l);

    int tp_size;
    fgets(temp_l,254,fp1);
    sscanf(temp_l,"tp_size=%d",&tp_size);
    int tpi_size;
    fgets(temp_l,254,fp1);
    sscanf(temp_l,"tpi_size=%d",&tpi_size);

    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"tp=");
    int* tp=read_int_vec(token,tp_size);
    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"tpi=");
    int* tpi=read_int_vec(token,tpi_size);
    
    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"c=");
    double* c=read_double_vec(token,l);
    for (int i = 0; i < l; i++) {
        c[i] = c[i]*60*dt/1000;
    }
    
    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"gamma=");
    double* b_gamma=read_double_vec(token,l);
    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"rho=");
    double* rho=read_double_vec(token,l);
    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"psi=");
    double* psi=read_double_vec(token,l*l);
    
    double eps_min;
    fgets(temp_l,254,fp1);
    sscanf(temp_l,"eps_min=10^%lf",&eps_min);
    eps_min=pow(10,eps_min);

    fclose(fp1);

    for (int i = 0; i < l*l; i++) {
        psi[i] = psi[i]/dt;
    }
    double psi1[l][l];
    
    for (int i=0;i<l;i++){
        for (int j=0;j<l;j++){
            psi1[i][j]=psi[i+j*l];
        }
    }
    temp=0;
    for (int i = 0; i < l; i++) {
        temp =temp + c[i];
    }
    double r[N];
    for (int i = 0; i < N; i++) {
        r[i] = temp;
    }

    fp1 = fopen("CBEC_2d_req.txt", "r");

    int K;
    fgets(temp_l,254,fp1);
    sscanf(temp_l,"K=%d",&K);
    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"Ki=");
    int* Ki=read_int_vec(token,K);
    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"q=");
    double* q=read_double_vec(token,K);

    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"s=");
    int* s=read_int_vec(token,K);
    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"d=");
    int* d=read_int_vec(token,K);
    fgets(temp_l,254,fp1);
    token=strtok(temp_l,"eps=");
    double* eps=read_double_vec(token,K);

    double alfa[K];
    double beta[K];
    for (int i = 0; i < K; i++) {
        alfa[i]=1;
        beta[i]=1;
    }
    double j_funct[4]={1,1000,10,1};
    for (int i = 0; i < K; i++) {
        q[i] = q[i]*60*dt/1000;
        d[i] = ceil(d[i]*60/dt);
        //eps[i] = eps[i]*0.8;
    }
    
    double Dt[K];
    double Dv[K];
    for (int i = 0; i < K; i++) {
        Dt[i]=N-s[i]-eps[i]*d[i]/dt;
        Dv[i]=(1-eps[i])*q[i]*d[i];
        if (Dt[i]<s[i]-1){
            Dt[i]=s[i]-1;
        }
        if (Dt[i]<1){
            Dt[i]=1;
        }
    }
    
    
    
    int delay1[l];
    int Ki1[K];
    int Ii1[l];
    for (int i=0;i<l;i++){
        delay1[i]= (int)delay[i];
        Ii1[i]=(int)Ii[i];
    }
    
    for (int i=0;i<K;i++){
        Ki1[i]= (int)Ki[i];
    }
    fclose(fp1);

    
    
    int NN=4*l*N+l*l*N*M+l*N*M+2*K*N+2*K;

    double *x0= initial_solution( N,  M,  K,  l, tp_size, tpi_size, tp, tpi, c, b_gamma, rho, psi1, delay1, Ki1, Ii1, q, s, d,  eps, eps_min,total_delay,H0,P0,L0);
    
    fp1 = fopen("temp.txt", "w");
    for (int i=0;i<NN;i++){
        fprintf(fp1,"%f\n",x0[i]);
    }
    fclose(fp1);

    
    emxArray_real_T *f;
    emxArray_real_T *A_beg;
    emxArray_real_T *A_nnz;
    emxArray_int32_T *A_indx;
    emxArray_real_T *A_val;
    emxArray_char_T *A_sense;
    emxArray_real_T *b;
    emxArray_real_T *Aeq_i;
    emxArray_real_T *Aeq_j;
    emxArray_real_T *Aeq_d;
    emxArray_real_T *beq;
    emxArray_real_T *lb;
    emxArray_real_T *ub;
    emxArray_char_T *ctype;
    emxArray_real_T *A_j;
    
    emxInitArray_real_T(&f, 1);
    emxInitArray_real_T(&A_beg, 1);
    emxInitArray_real_T(&A_nnz, 1);
    
    emxInitArray_int32_T(&A_indx, 1);
    //emxInitArray_real_T(&A_indx, 1);
    emxInitArray_real_T(&A_val, 1);
    emxInitArray_char_T(&A_sense, 2);
    emxInitArray_real_T(&b,1);
    emxInitArray_real_T(&Aeq_i, 1);
    emxInitArray_real_T(&Aeq_j, 1);
    emxInitArray_real_T(&Aeq_d, 1);
    emxInitArray_real_T(&beq, 1);
    emxInitArray_real_T(&lb, 1);
    emxInitArray_real_T(&ub, 1);
    emxInitArray_char_T(&ctype, 2);
    emxInitArray_real_T(&A_j, 1);
    
    double A_r;
    
    
     get_MILP(N,M,K,l,tp_size,tpi_size,total_delay,P0,L0,H0,tp,tpi,c,r,b_gamma,rho,psi,delay,eps_min,Ki,Ii,q,s,d,eps,alfa,beta,Dt,Dv,j_funct,
             f,A_beg,A_nnz,A_indx,A_val,A_sense,b,lb,ub,ctype,A_j,&A_r);
    
    int nnz_ele=A_val->size[0];    int numcols=f->size[0]; int numrows=b->size[0];
    double  objective[numcols];
    double  lp_lb[numcols];
    double  lp_ub[numcols];
    double  rhs[numrows];
    char  sense[numrows];
    int  matbeg[numcols];
    int matcnt[numcols];
    //int matind[nnz_ele];
    //double*  matval=(double*) malloc(nnz_ele*sizeof(double*));
    //matval[0]=1;
    for(int i=0;i<numcols;i++){
        objective[i]=f->data[i];
        lp_lb[i]=lb->data[i];
        lp_ub[i]=ub->data[i];
        matbeg[i]=(int)A_beg->data[i];
        matcnt[i]=(int)A_nnz->data[i];
        //printf("%d\t%f\n",i,lp_ub[i]);
        
    }
    
    for(int i=0;i<numrows;i++){
        rhs[i]=b->data[i];
        sense[i]=A_sense->data[i];
    }
    
    
    //Create the CPLEX problem
    int status = 0;
    int* beg=(int*)malloc(1 * sizeof(int*));
    beg[0]=0;
    int* varindices=(int*)malloc(NN * sizeof(int*));;
    for (int i=0;i<NN;i++){
        varindices[i]=i;
    }
    CPXENVptr env = CPXopenCPLEX (&status);
    CPXLPptr lp= CPXcreateprob (env, &status, "myprob");
    status= CPXcopylp (env, lp, numcols, numrows, 1, objective, rhs,
                       sense, matbeg, matcnt, A_indx->data, A_val->data,
                       lp_lb, lp_ub, NULL);
    status = CPXcopyctype (env, lp, ctype->data);
    
    status= CPXaddmipstarts(env, lp, 1, NN, beg, varindices,
                            x0, NULL, NULL);
    status = CPXsetintparam (env, CPXPARAM_ClockType, 1);
    status = CPXsetdblparam (env, CPXPARAM_TimeLimit, 10.0);
    
    status = CPXsetintparam (env, CPXPARAM_ScreenOutput, CPX_ON);
    status = CPXwriteprob(env, lp,"my_probv1.lp",NULL);

    status = CPXmipopt (env, lp);
    
    
    
/*    const int N=36;
    int M=20;
    const int l=7;
    int dt=60;
    const double V0[l]={0,0,0,0,0,0,0};
    const double R0[l]={0,0,0,0,0,0,0};
    const double H0[l]={0,0,0,0,0,0,0};
    //const double b_gamma[l]={0.1,0.1,0.1,0.1,0.1,0.1,0.1};
    const double b_gamma[l]={0.03,0.03,0.03,0.03,0.03,0.03,0.03};
    const double rho[l]={0.4,0.4,0.4,0.4,0.4,0.4,0.4};

    const int tp_size=15;
    const int tpi_size=12;
    int tp[tp_size]={8,9,10,11,12,20,21,22,23,24,32,33,34,35,36};
    int tpi[tpi_size]={9,10,11,12,21,22,23,24,33,34,35,36};

    for (int i=0;i<tp_size;i++){
        tp[i]=tp[i]-1;
    }
    for (int i=0;i<tpi_size;i++){
        tpi[i]=tpi[i]-1;
    }

    double c[l]={320,200,110,75,240,75,100};
    for (int i = 0; i < l; i++) {
        c[i] = c[i]*60*dt/1000;
    }
    
    double temp=0;
    for (int i = 0; i < l; i++) {
        temp =temp + c[i];
    }
    double r[N];
    for (int i = 0; i < N; i++) {
        r[i] = temp;
    }
    double psi[l*l]= {45,45,60,60,75,90,90,45,30,45,45,60,75,75,60,45,30,45,60,75,75,60,45,45,30,60,75,75,75,60,60,60,30,60,60,90,75,75,75,60,20,20,90,75,75,75,60,30,20};
    
    for (int i = 0; i < l*l; i++) {
            psi[i] = psi[i]/dt;
        }
    double psi1[l][l];
    
    for (int i=0;i<l;i++){
        for (int j=0;j<l;j++){
            psi1[i][j]=psi[i+j*l];
        }
    }

    //const double tau[l]={120,30,60,90,60,120,120};
    const double tau[l]={180,60,60,60,60,60,60};
    double delay[l];
    for (int i = 0; i < l; i++) {
        temp=ceil(tau[i]/dt);
        if (temp<1){
            delay[i]=1;
        }
        else{
            delay[i]=temp;
        }
    }
    const double eps_min=pow(10,-8);


    const int K=33;
    const double Ki[K]={1,7,2,3,6,4,5,1,4,1,2,3,6,4,6,7,4,4,7,4,4,4,1,6,4,4,7,4,4,4,6,6,4};
    const double Ii[l]={0,1,1,3,3,5,5};
    double q[K]={20,20,20,20,13,13.75,25,20,25,20,20,20,13,25,15,30,12,12,30,25,25,25,20,13,25,25,30,12,25,25,15,15,25};
    double    s[K]={1,1,1,1,1,1,1,6,13,13,13,13,13,13,13,13,13,13,17,18,18,25,25,25,25,25,25,25,25,28,29,29,30};
    double    d[K]={8,8,8,8,8,8,2,6,8,8,8,4,8,8,7,8,8,8,8,4,6,8,8,6,6,4,8,4,7,8,8,8,5};
    double eps[K]= {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    double alfa[K]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    double beta[K]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    double j_funct[4]={1,1000,10,1};
    for (int i = 0; i < K; i++) {
        q[i] = q[i]*60*dt/1000;
        d[i] = ceil(d[i]*60/dt);
        eps[i] = eps[i]*0.8;
        //s[i]=s[i]*4; //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    }
    
    double Dt[K];
    double Dv[K];
    for (int i = 0; i < K; i++) {
        Dt[i]=N-s[i]-eps[i]*d[i]/dt;
        Dv[i]=(1-eps[i])*q[i]*d[i];
        if (Dt[i]<s[i]-1){
            Dt[i]=s[i]-1;
        }
        if (Dt[i]<1){
            Dt[i]=1;
        }
    }
    

    
    int delay1[l];
    int Ki1[K];
    int Ii1[l];
    for (int i=0;i<l;i++){
        delay1[i]= (int)delay[i];
        Ii1[i]=(int)Ii[i];
    }

    for (int i=0;i<K;i++){
        Ki1[i]= (int)Ki[i];
    }
    
    
    double *x0= initial_solution( N,  M,  K,  l, tp_size, tpi_size, tp, tpi, c, b_gamma, rho, psi1, delay1, Ki1, Ii1, q, s, d,  eps, eps_min);
    int NN=4*l*N+l*l*N*M+l*N*M+2*K*N+2*K;
    FILE* fp1;
    fp1 = fopen("temp.txt", "w");
    for (int i=0;i<NN;i++){
        fprintf(fp1,"%f\n",x0[i]);
    }
    fclose(fp1);

    emxArray_real_T *f;
    emxArray_real_T *A_beg;
    emxArray_real_T *A_nnz;
    emxArray_int32_T *A_indx;
    emxArray_real_T *A_val;
    emxArray_char_T *A_sense;
    emxArray_real_T *b;
    emxArray_real_T *Aeq_i;
    emxArray_real_T *Aeq_j;
    emxArray_real_T *Aeq_d;
    emxArray_real_T *beq;
    emxArray_real_T *lb;
    emxArray_real_T *ub;
    emxArray_char_T *ctype;

    emxInitArray_real_T(&f, 1);
    emxInitArray_real_T(&A_beg, 1);
    emxInitArray_real_T(&A_nnz, 1);
    
    emxInitArray_int32_T(&A_indx, 1);
    //emxInitArray_real_T(&A_indx, 1);
    emxInitArray_real_T(&A_val, 1);
    emxInitArray_char_T(&A_sense, 2);
    emxInitArray_real_T(&b, 1);
    emxInitArray_real_T(&Aeq_i, 1);
    emxInitArray_real_T(&Aeq_j, 1);
    emxInitArray_real_T(&Aeq_d, 1);
    emxInitArray_real_T(&beq, 1);
    emxInitArray_real_T(&lb, 1);
    emxInitArray_real_T(&ub, 1);
    emxInitArray_char_T(&ctype, 2);

    
    get_MILP(N,M,K,l,dt,tp_size,tpi_size,V0,R0,H0,tp,tpi,c,r,b_gamma,rho,psi,delay,eps_min,Ki,Ii,q,s,d,eps,alfa,beta,Dt,Dv,j_funct,
             f,A_beg,A_nnz,A_indx,A_val,A_sense,b,lb,ub,ctype);
    
//    double *Alp_beg=A_beg->data;
//    double *Alp_nnz=A_nnz->data;
//    double *Alp_indx=A_indx->data;
    int nnz_ele=A_val->size[0];    int numcols=f->size[0]; int numrows=b->size[0];
    double  objective[numcols];
    double  lp_lb[numcols];
    double  lp_ub[numcols];
    double  rhs[numrows];
    char  sense[numrows];
    int  matbeg[numcols];
    int matcnt[numcols];
    //int matind[nnz_ele];
    //double*  matval=(double*) malloc(nnz_ele*sizeof(double*));
    //matval[0]=1;
    for(int i=0;i<numcols;i++){
        objective[i]=f->data[i];
        lp_lb[i]=lb->data[i];
        lp_ub[i]=ub->data[i];
        matbeg[i]=(int)A_beg->data[i];
        matcnt[i]=(int)A_nnz->data[i];
        //printf("%d\t%f\n",i,lp_ub[i]);

    }

    for(int i=0;i<numrows;i++){
        rhs[i]=b->data[i];
        sense[i]=A_sense->data[i];
    }

    
    //Create the CPLEX problem
    int status = 0;
    int* beg=(int*)malloc(1 * sizeof(int*));
    beg[0]=0;
    int* varindices=(int*)malloc(NN * sizeof(int*));;
    for (int i=0;i<NN;i++){
        varindices[i]=i;
    }
    CPXENVptr env = CPXopenCPLEX (&status);
    CPXLPptr lp= CPXcreateprob (env, &status, "myprob");
    //status= CPXcopylp (env, lp, numcols, numrows, 1, objective, rhs,
    //                  sense, matbeg, matcnt, matind, matval,
    //                  lp_lb, lp_ub, NULL);
    status= CPXcopylp (env, lp, numcols, numrows, 1, objective, rhs,
                                         sense, matbeg, matcnt, A_indx->data, A_val->data,
                                         lp_lb, lp_ub, NULL);
    status = CPXcopyctype (env, lp, ctype->data);

    status= CPXaddmipstarts(env, lp, 1, NN, beg, varindices,
                             x0, NULL, NULL);
    status = CPXsetintparam (env, CPXPARAM_ClockType, 1);
    status = CPXsetdblparam (env, CPXPARAM_TimeLimit, 10.0);
    
    status = CPXsetintparam (env, CPXPARAM_ScreenOutput, CPX_ON);
//    status = CPXmipopt (env, lp);
    //status = CPXlpopt(env, lp);
    status = CPXwriteprob(env, lp,"my_probv1.lp",NULL);

    //printf('%f',n_c[0]);
    printf("Hello, World!\n");
    return 0;
 */
}



