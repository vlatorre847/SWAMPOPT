//
//  my_functions.c
//  initial_solution
//
//  Created by Valyria on 25/09/20.
//  Copyright © 2020 Valyria. All rights reserved.
//

#include "my_functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "my_math.h"
#include <math.h>
#include <string.h>


int* read_int_vec(char* token, int v_length){
    int* vec=(int*)malloc(v_length*sizeof(int*));
    int index_s=0;
    char temp_token;
    int index_e=0;
    for (int i=0;i<v_length;i++){
        temp_token=token[index_s];
        while (temp_token!=' '){
            index_e=index_e+1;
            temp_token=token[index_e];
        }
        char token1[254]="";
        for (int j=index_s;j<index_e;j++){
            token1[j-index_s]=token[j];
        }
        sscanf(token1,"%d",&vec[i]);
        //printf("%d\n",vec[i]);
        index_s=index_e+1;
    }
    return vec;
}

double* read_double_vec(char* token, int v_length){
    double* vec=(double*)malloc(v_length*sizeof(double*));
    int index_s=0;
    char temp_token;
    int index_e=0;
    for (int i=0;i<v_length;i++){
        temp_token=token[index_s];
        while (temp_token!=' '){
            index_e=index_e+1;
            temp_token=token[index_e];
        }
        char token1[254]="";
        for (int j=index_s;j<index_e;j++){
            token1[j-index_s]=token[j];
        }
        sscanf(token1,"%lf",&vec[i]);
        //printf("%lf\n",vec[i]);
        index_s=index_e+1;
    }
    return vec;
}

void check_mat(int n_rows, int n_cols, double mat[n_rows][n_cols]){
    for(int i=0;i<n_rows;i++){
        for (int j=0;j<n_cols;j++){
            printf("%f\t",mat[i][j]);
        }
        printf("\n");
    }
}

void check_p(int n_rows, int n_cols, double **mat){
    for(int i=0;i<n_rows;i++){
        for (int j=0;j<n_cols;j++){
            //printf("%f\t",mat[i+n_rows*j]);
            printf("%f\t",(double)mat[i][j]);
        }
        printf("\n");
    }
}


void check_int(int n_rows, int n_cols, int **mat){
    for(int i=0;i<n_rows;i++){
        for (int j=0;j<n_cols;j++){
            //printf("%f\t",mat[i+n_rows*j]);
            printf("%d\t",mat[i][j]);
        }
        printf("\n");
    }
}


bool my_ismember(int n,double a,int arr[n]){
    bool res=false;
    for(int i= 0; i<n;i++){
        if (arr[i]==a){
            res=true;
            break;
        }
    }
    return res;
}

int* my_find(int n,int val, int arr[n], int ops, int* size_arg){
    /*ops=-1, take the first element, ops=1, take the last element*/
    int* result=(int*)malloc(n * sizeof(int*));
    int size=0;
    for( int i=0;i<n;i++){
        if (arr[i]==val){
            result[size]=i;
            size=size+1;
        }
    }
    int* result1=(int*)malloc( sizeof(int*));
    if (ops==-1){
        result1=&result[0];
        if (size>0){
        size=1;
        }
        *size_arg=size;
        return result1;
    }
        else{
            if (ops==1){
            result1=&result[size-1];
            if (size>0){
                size=1;
            }
           *size_arg=size;
            return result1;
            }
        }
    *size_arg=size;
    return result;
}


int* find_greater(int n,int val, int arr[n], int ops, int* size_arg){
    /*ops=-1, take the first element, ops=1, take the last element*/
    int* result=(int*)malloc(n * sizeof(int*));
    int size=0;
    for( int i=0;i<n;i++){
        if (arr[i]>=val){
            result[size]=i;
            size=size+1;
        }
    }
    int* result1=(int*)malloc( sizeof(int*));
    if (ops==-1){
        result1=&result[0];
        size=1;
        *size_arg=size;
        return result1;
    }
    else{
        if (ops==1){
            result1=&result[size-1];
            size=1;
            *size_arg=size;
            return result1;
        }
    }
    *size_arg=size;
    return result;
}




double** chan_req(int N,int l,int K,int tp_size,double gamma[l],double q[K],int Ii[l],int Ki[K],int tp[tp_size],int delay[l],int** D){
    double** CR = (double**)malloc(l * sizeof(double**));
    
    for (int i = 0; i < l; i++) {
        CR[i] = (double*)malloc((N+1) * sizeof(double*));
    }
    
    for (int i=0;i<l;i++){
        int* ir_index=(int*)malloc(K * sizeof(int*));
        int ir_index_size=0;
        for( int ii=0;ii<K;ii++){
            if (Ki[ii]==i+1){
                ir_index[ir_index_size]=ii;
                ir_index_size=ir_index_size+1;
            }
        }
        for(int j=delay[i];j<N;j++){
            if (my_ismember(tp_size, j-delay[i]+1,tp)){
                CR[i][j-delay[i]]=CR[i][j-delay[i]-1];
            }else{
                double qd[ir_index_size];
                for (int k=0;k<ir_index_size;k++){
                    qd[k]=q[ir_index[k]]*D[ir_index[k]][(j-1+delay[i])];
            }
                CR[i][j-delay[i]]=my_sum(ir_index_size,qd)/pow(1-gamma[i], delay[i]);
            }
        }
    }
    int prec[l];
    for (int ij=0;ij<l;ij++){
        prec[ij]=Ii[ij]-1;
    }

    for (int i=0;i<l;i++){
        int temp=delay[i]-1;
        int prec_c=prec[i];
        while (prec_c!=-1){
            temp=temp+delay[prec_c];
            prec_c=prec[prec_c];
        }
        for(int j=temp;j<N;j++){
            if (CR[i][j]==0)
                continue;
            prec_c=prec[i];
            int jj=0;
            double appo=CR[i][j];
            while (prec_c!=-1){
                jj=jj+delay[prec_c];
                appo=appo/pow(1-gamma[prec_c],jj);
                CR[prec_c][j-jj]=CR[prec_c][j-jj]+appo;
                prec_c=prec[prec_c];
            }
        }
    }
    return CR;
}

double gk_workload(int l,int N,int j,int ii,double psi[l][l],int** G){
    double psi_temp=0;
    double psi_temp1=0;
    if (ii==1&&j==28){
        
    }
    int jj=0;
    int G_temp[l];
    for (int i=0;i<l;i++){
        G_temp[i]=(int)G[i][j-jj];
    }
    int size_arg_last=0;
    int last_gate=*my_find( l, 1, G_temp,1,&size_arg_last);
    while (size_arg_last==0){
        jj=jj+1;
        for (int i=0;i<l;i++){
            G_temp[i]=(int)G[i][jj];
        }
        last_gate=*my_find( l, 1, G_temp,1,&size_arg_last);
        if (jj>j)
            return 0;
    }
    
    //Until here founds the last gate that has been activated
    
    int G_appo[l];
    for (int i=0;i<l;i++){
        G_appo[i]=(int)G[i][j];
    }
    int size_arg=0;
    int* index1=my_find( l, 1, G_appo,0,&size_arg);
    
    if (size_arg==0){
        if (psi[last_gate][ii]<=1){
            psi_temp=psi[last_gate][ii];
        }else
            psi_temp=-jj+1+psi[last_gate][ii];
    }else{
        G_appo[ii]=1;
        index1=my_find( l, 1, G_appo,0,&size_arg);
        if (psi[last_gate][index1[0]]<=1){
            psi_temp=psi[last_gate][index1[0]];
        }else
            psi_temp=0;
        for(int iii=1;iii<size_arg;iii++){
            psi_temp=psi_temp+psi[index1[iii-1]][index1[iii]];
        }
    }
    int jn=j;
    while (jn+1!=N){
        int G_appo1[l];
        for (int i=0;i<l;i++){
            //G_appo1[i]=G[i+(jn+1)*l];
            G_appo1[i]=G[i][jn+1];
        }
        int sum_g=sum_i(l,G_appo1);
        if(sum_g==0){
            jn=jn+1;
            continue;
        }
        int* index2=my_find(l, 1, G_appo1, 0, &size_arg);
        if (psi[ii][index2[0]]<=1){
                psi_temp1=psi[ii][index2[0]];
            for (int iii=1;iii<size_arg;iii++){
                psi_temp1=psi_temp1+psi[index2[iii-1]][index2[iii]];
            }
        }else
            psi_temp1=j-jn+psi[ii][index2[0]];
        psi_temp=my_max(psi_temp,psi_temp1);
        break;
    }
    return psi_temp;
}

void adjust_flow(int N,int M,int l,int ii,int j,int tp_size,double** P,int** G,double** H,double** CR,double c[l],double psi[l][l],int tp[tp_size],double rho[l], double psi_temp, int total_delay, double V0[total_delay],int delay[l]){
    int sum_delay=0;
    double tempP=0;
    if (j<delay[ii]){
        for (int iii=0;iii<delay[ii];iii++){
            tempP=V0[sum_delay+iii];
            sum_delay=sum_delay+delay[ii];
        }
    }else{
        tempP=P[ii][j-1];
    }

    if (P[ii][j]<rho[ii]*H[ii][j]*c[ii]){
        if (my_ismember(tp_size, j, tp)||psi_temp>1){
            P[ii][j]=tempP;//P[ii][(j-1)];
            if (P[ii][j-1]>CR[ii][j])
                CR[ii][j]=P[ii][j-1];
        }else{
            
            G[ii][j]=1;
            int jj=j+1;
            while (jj<N){
                double psi_temp1=gk_workload(l, N, jj, ii, psi, G);
                if (psi_temp1<=1 || G[ii][jj]==1){
                    if (my_ismember(tp_size, jj, tp)){
                        jj=jj+1;
                        continue;
                    }
                    
                    G[ii][jj]=1;
                    break;
                }
                jj=jj+1;
            }
            double minp=0;
            for (int jjj=j;jjj<=jj-1;jjj++){
                minp=my_min(minp,P[ii][jjj]/(c[ii]*rho[ii]));
            }
            for (int jjj=j;jjj<=jj-1;jjj++){
                H[ii][jjj]=minp;
            }
        }
    }else{
        
        
        if (P[ii][j]>H[ii][j]*c[ii]){
            if (my_ismember(tp_size, j, tp)||psi_temp>1){
                P[ii][j]=tempP;//P[ii][(j-1)];
                if (P[ii][(j-1)]>CR[ii][j])
                    CR[ii][j]=P[ii][j-1];
            }else{
                
                G[ii][j]=1;
                int jj=j+1;
                while (jj<N){
                    double psi_temp1=gk_workload(l, N, jj, ii, psi, G);
                    //printf("%f\n",G[ii][jj]);
                    if (psi_temp1<=1 || G[ii][jj]==1){
                        if (my_ismember(tp_size, jj, tp)){
                            jj=jj+1;
                            continue;
                        }
                        if(ii==3&&jj==25){
                            //printf("qui");
                        }
                        
                        G[ii][jj]=1;
                        break;
                    }
                    jj=jj+1;
                }
                double maxp=0;
                for (int jjj=j;jjj<=jj-1;jjj++){
                    maxp=my_max(maxp,P[ii][jjj]/c[ii]);
                }
                for (int jjj=j;jjj<=jj-1;jjj++){
                    H[ii][jjj]=maxp;
                }
            }
        }
    }

}


void add_irrigation(int N,int M,int l,int K,int tp_size,double** L,double** P,int** G,double** H,int** D,int stopj[l],int startj[l],double q[K],double c[l],int Ii[l],int Ki[K],double gamma[l],double psi[l][l],int tp[tp_size],double rho[l],int delay[l], const double eps_min,int total_delay,double V0[total_delay],double R0[l]){
    
    for (int j=startj[0];j<N;j++){
        P[0][j]=c[0];
    }
    //check_p(l,N,L);
    
    double** CR=chan_req( N, l, K, tp_size, gamma, q, Ii, Ki, tp, delay, D);
    //check_p(l,N,G);
    for (int i=0;i<l;i++){
        int size_ir=0;
        int* ir_index=my_find(K, i+1, Ki, 0, &size_ir);
        int size_temp_ch=0;
        int* temp_ch=my_find(l, i+1, Ii, 0, &size_temp_ch);
        double qd[size_ir];
        double tempL=0;
        double tempP=0;

        for (int j=0;j<N;j++){
            if (j==0){
                tempL=R0[i];
            }else{
                tempL=L[i][j-1];
            }
            int sum_delay=0;
            if (j<delay[i]){
                for (int ii=0;ii<delay[i];ii++){
                    tempP=V0[sum_delay+ii];
                    sum_delay=sum_delay+delay[i];
                }
            }else{
                tempP=P[i][j-delay[i]];
            }
            
            
            
            for (int kk=0;kk<size_temp_ch;kk++){
                int ii=temp_ch[kk];
                double psi_temp= gk_workload(l, N, j, ii, psi, G);
                for (int k=0;k<size_ir;k++){
                    qd[k]=q[ir_index[k]]*D[ir_index[k]][j];
                }
                double sum_cr=0;
                for (int k=0;k<size_temp_ch;k++){
                    if (temp_ch[k]!=ii)
                        sum_cr=sum_cr+CR[temp_ch[k]][j];
                }
                double tempf=pow(1-gamma[i],delay[i])*tempP+(1-gamma[i])*tempL-my_sum(size_ir,qd)-sum_cr;
                int temp_size;
                int* down_ch=my_find(l, ii+1, Ii, 0, &temp_size);
                if (temp_size!=0)
                    P[ii][j]=my_min(tempf,c[ii]);
                else{
                    int size_ir1;
                    int* ir_index1=my_find(K, ii+1, Ki, 0, &size_ir1);
                    double qd1[size_ir1];
                    for (int k=0;k<size_ir1;k++){
                        qd1[k]=q[ir_index1[k]]*D[ir_index1[k]][(int)my_min(j+delay[ii],N)];
                    }
                    double tempd=my_sum(size_ir1,qd1)/pow(1-gamma[ii],delay[ii]);
                    P[ii][j]=my_min(tempf,my_min(c[ii],tempd));
   //                 if (j>=stopj[ii]||j<startj[ii])
     //                   P[ii][j]=0;
                }
                
/**********************************************************************************
 %Set the values for H
**********************************************************************************/
                adjust_flow( N, M, l, ii, j, tp_size, P, G, H, CR, c, psi, tp,rho, psi_temp, total_delay,V0,delay);
            }
            /**********************************************************************************
             %Set the values for L
             **********************************************************************************/
            
            double sum_p=0;
            for (int k=0;k<size_temp_ch;k++){
                sum_p=sum_p+P[temp_ch[k]][j];
            }
            double qd1[size_ir];
            for (int k=0;k<size_ir;k++){
                qd1[k]=q[ir_index[k]]*D[ir_index[k]][j];
            }

            L[i][j]=pow(1-gamma[i],delay[i])*tempP+(1-gamma[i])*tempL-my_sum(size_ir,qd1)-sum_p;
            if (L[i][j]<eps_min&&-L[i][j]<eps_min)
                L[i][j]=0;

        }
        for (int ii=1;ii<l;ii++){
            for (int jj=1;jj<N;jj++){
                if (H[ii][jj]==H[ii][(jj-1)]&&G[ii][jj]==1){
                    G[ii][jj]=0;
                }
            }
        }
    }
}
