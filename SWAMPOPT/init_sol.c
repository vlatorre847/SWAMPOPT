//
//  init_sol.c
//  initial_solution
//
//  Created by Valyria on 25/09/20.
//  Copyright © 2020 Valyria. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include "init_sol.h"
#include "my_math.h"
#include "my_functions.h"

double*  initial_solution(const int N, const int M, const int K, const int l,const int tp_size,const int tpi_size, const int                                 tp[tp_size], const int tpi[tpi_size], const double c[l], const double b_gamma[l], const double rho[l], const double psi[l][l],const int delay[l], const int Ki[K],const int Ii[l], const double q[K], const int s[K], const int d[K], const double eps[K], const double eps_min, int  total_delay,const double H0[l], const double V0[total_delay], const double R0[l])
{
    /*
    N: Number of time intervals
    M: Number of operations
    K: Numer of irrigations
    l: Number of channels
    Ki: Set of the sets of off-takes on the channels
    Ii: Set of the sets of the channels downstream every channel
    tp_size: size of tp
    tpi_size: size of tpi
    tp: Time intervals where the gate-keeper cannot operate
    tpi: Time intervals where irrigations cannot start
    q: Quantity of water required by the off-take for time interval
    s: Desidered starting time interval for the irrigation
    d: Desidered duration for the irrigation
    eps: minimum ration of water that must be delivered for irrigation
    b_gamma: seepage ration for channel
    psi: travelling time from one gate to another
    c: Maximum inlet volume capacity for every gate
    rho: Minimum ratio of water that must flow through an open gate
     */
 
// Variables Declaration
    double** P= (double**)malloc(l * sizeof(double**));
    double** L= (double**)malloc(l * sizeof(double**));
    int** G= (int**)malloc(l * sizeof(int**));
    double** H= (double**)malloc(l * sizeof(double**));
    int** D= (int**)malloc(K * sizeof(int**));
     
    double F[l][l][M][N];
    memset(F,0,sizeof(double)*l*l*M*N);
    double E[l][M][N];
    memset(E,0,sizeof(double)*l*M*N);
    double S[K][N];
    memset(S,0,sizeof(double)*K*N);

    for (int i = 0; i < l; i++) {
        P[i] = (double*)malloc((N) * sizeof(double*));
        L[i] = (double*)malloc((N) * sizeof(double*));
        G[i] = (int*)malloc((N) * sizeof(int*));
        H[i] = (double*)malloc((N) * sizeof(double*));
    }
    for (int i = 0; i < K; i++) {
        D[i] = (int*)malloc((N) * sizeof(int*));
    }
    
    // Variables allocation
    for (int i = 0; i < l; i++) {
        for (int j = 0; j < N; j++) {
            P[i][j]=0;
            L[i][j]=0;
            G[i][j]=0;
            H[i][j]=0;
        }
    }
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < N; j++) {
            D[i][j]=0;
        }
    }

    int dmin[K]; //Minimum duration for irrigation
    for(int i=0;i<K;i++){
        dmin[i]=ceil(d[i]*eps[i]);
    }
    int startj[l];
    /*%%%%%%%%P,L,S,D in the first (source) channel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    
    startj[0]=ceil(psi[0][0])-1;
    G[0][startj[0]]=1;
    for (int i=0; i<N; i++) {
        P[0][i]=c[0];
    }
    int* ir_index=(int*)malloc(K * sizeof(int*));;
    int ir_index_size=0;
    
    for( int i=0;i<K;i++){
        if (Ki[i]==1){//In Ki the first channel is the number 1, that has index 0
            ir_index[ir_index_size]=i;
            ir_index_size=ir_index_size+1;
        }
    }
    
    for(int ik=0; ik<ir_index_size;ik++){
        int j= my_max(delay[0]+startj[0],s[ir_index[ik]]-1);
        if (j+dmin[ir_index[ik]]>N)
            continue;
        while (j<N){
            double qD[ir_index_size];
            for (int k=0;k<ir_index_size;k++){
                qD[k]=q[ir_index[k]]*D[ir_index[k]][j-1+delay[0]];
            }
            double sum_qD=my_sum(ir_index_size,qD);
            if (P[0][j-1]+L[0][j-1]-sum_qD<q[ir_index[ik]]){
                j=j+1;
                continue;
            }
            S[ir_index[ik]][j]=1;
            for (int ii=j;ii<j+dmin[ir_index[ik]];ii++){
                D[ir_index[ik]][ii]=1;
            }
            break;
        }
    }
    int stopj[l];
    for (int i=0;i<l;i++){
        stopj[i]=N+1;
    }
    /*For cycle for opening the channels, so that the time necessary to do the
     operations does not exceede the duration of the time interval itself
     The gates are opened following their sequence
     (gate 1 is opened first, then gate 2 and so on)*/
    int j_temp=1;
    double temppsi=0;
    for (int i=1;i<l;i++){
        while (j_temp<N-1){
            if (temppsi+psi[i-1][i]>1){
                j_temp=j_temp+1;
                temppsi=temppsi-1;
                continue;
            }
            startj[i]=j_temp;
            if (psi[i-1][i]>1)
                temppsi=0;
                else
                temppsi=temppsi+psi[i-1][i];
            break;
        }
    }
    
    //for (int jj=0;jj<l;jj++){
    //    printf("%d\n",startj[jj]);
    //}
    
    for (int j=0;j<N;j++){
        H[0][j]=P[0][j]/c[0];
    }
    int sum_delay=0;
    for (int i=1;i<l;i++){
        for (int ii=0;ii<delay[i];ii++){
            P[i][ii]=V0[sum_delay+ii];
            sum_delay=sum_delay+delay[i];
        }
        for (int j=0;j<N;j++){
            H[i][j]=H0[i];
        }
    }
    //check_p(l,N,P);
    /*%%%%%%%S, D, G%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    add_irrigation(N, M, l, K, tp_size, L, P, G, H, D, stopj, startj, q, c, Ii, Ki, b_gamma, psi, tp, rho, delay, eps_min,total_delay,V0,R0);
    //check_p(l,N,P);
    int** dappo= (int**)malloc(l * sizeof(int**));
    for (int i = 0; i < l; i++) {
        dappo[i] = (int*)malloc((N) * sizeof(int*));
    }

    for (int ii=1;ii<l;ii++){
        int size_ir;
        int* ir_index=my_find(K, ii+1, Ki, 0, &size_ir);
        int  dmin_ir[size_ir];
        for (int jj=0;jj<size_ir;jj++){
            dmin_ir[jj]=dmin[ir_index[jj]];
        }
        int* indexk= my_sort(dmin_ir, size_ir,1);
        for (int jj=0;jj<size_ir;jj++){
            //indexk[jj]=ir_index[indexk[jj]];
        }
        for (int ik=0;ik<size_ir;ik++){
            bool flag=false;
            int t=indexk[ik];
            int j=my_max(startj[ii]+1,s[ir_index[t]]);
            while (j<N){
                if (j+dmin[ir_index[t]]>N-1&&flag==0){
                    flag=1;
                    j=startj[ii]+1;
                }
                if (j+dmin[ir_index[t]]>N && flag==1) //If you cannot do the irrigation in the considered time interval break
                    break;
                if (my_ismember(tpi_size, j, tpi)){
                    j=j+1;
                    continue;
                }
                // Creation of the temporary variables to see if the irrigation is feasible
                double** P1= (double**)malloc(l * sizeof(double**));
                double** L1= (double**)malloc(l * sizeof(double**));
                int**    G1= (int**)malloc(l * sizeof(int**));
                double** H1= (double**)malloc(l * sizeof(double**));
                int**    D1=(int**)malloc(K * sizeof(int**));
                
                for (int i = 0; i < l; i++) {
                    P1[i] = (double*)malloc((N) * sizeof(double*));
                    L1[i] = (double*)malloc((N) * sizeof(double*));
                    G1[i] = (int*)malloc((N) * sizeof(int*));
                    H1[i] = (double*)malloc((N) * sizeof(double*));
                    for (int jj=0;jj<N;jj++){
                        P1[i][jj]=P[i][jj];
                        L1[i][jj]=L[i][jj];
                        G1[i][jj]=G[i][jj];
                        H1[i][jj]=H[i][jj];
                    }
                }
                
                for (int i = 0; i < K; i++) {
                    D1[i] = (int*)malloc((N) * sizeof(int*));
                    for (int jj=0;jj<N;jj++){
                    D1[i][jj]=D[i][jj];
                    }
                }
                for(int i = j; i < j+dmin[ir_index[t]]; i++) {
                    D1[ir_index[t]][i]=1;
                }
                if (ii==3&&ik==9){
                    //check_int(l,N,G);
                    //check_int(l,N,G1);
                    //printf("\n");
                }
                //if(ii==6){
                 //   printf("qui");
                //check_int(l,N,G);
                //}
                add_irrigation(N, M, l, K, tp_size, L1, P1, G1, H1, D1, stopj, startj, q, c, Ii, Ki, b_gamma, psi, tp, rho, delay, eps_min,total_delay,V0,R0);
                //check_p(l,N,P1);

                bool flag_ir=false;
                for (int i=0;i<l;i++){
                    for (int jj=0;jj<N;jj++){
                        if(L1[i][jj]<0||P1[i][jj]<0){
                            flag_ir=true;
                            break;
                        }
                    }
                }
                if (flag_ir){
                    j=j+1;
                    continue;
                }
                S[ir_index[t]][j]=1;

                for (int i = 0; i < l; i++) {
                    for (int jj=0;jj<N;jj++){
                        P[i][jj]=P1[i][jj];
                        L[i][jj]=L1[i][jj];
                        G[i][jj]=G1[i][jj];
                        H[i][jj]=H1[i][jj];
                    }
                }
                
                for (int i = 0; i < K; i++) {
                    for (int jj=0;jj<N;jj++){
                        D[i][jj]=D1[i][jj];
                    }
                }
//                if (G[5][7]==1)
//                    printf("%d,\n",G[5][7]);
                break;
            }
            for (int i=0;i<l;i++){
                int index1_size;
                int* ir_index1=my_find(K, i+1, Ki, 0, &index1_size);
                for (int jj=0;jj<N;jj++){
                    int D_temp=0;
                    for (int k=0;k<index1_size;k++){
                        D_temp=D_temp+D[ir_index1[k]][jj];
                    }
                    dappo[i][jj]=D_temp;
                }
            }
        }
        
        
        int dappo_t[N];
        for (int jj=0;jj<N;jj++){
            dappo_t[jj]=dappo[ii][jj];
        }
            
        int temp_size;
        int* temp=find_greater(N, 0, dappo_t, 1, &temp_size);
        stopj[ii]=*temp;
    }
        
    // %%%%%%%%E%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int im= nnz_int(l,N,G);
    int m=M-1;
    int k=1;
    int indexi[im];
    int indexm[im];
    int indexn[im];
    memset(indexi,0,sizeof(int)*im);
    memset(indexm,0,sizeof(int)*im);
    memset(indexn,0,sizeof(int)*im);
    //check_int(l,N,G);

    for (int jj=N-1;jj>=0;jj--){
        for (int i=l-1;i>=0;i--){
            if (G[i][jj]==1){
                E[i][m][jj]=1;
                indexi[im-k]=i;
                indexm[im-k]=m;
                indexn[im-k]=jj;
                k=k+1;
                m=m-1;
            }
        }
    }
    F[indexi[0]][indexi[0]][indexm[0]][indexn[0]]=1;
    for (int m=1;m<im;m++){
        F[indexi[m]][indexi[m-1]][indexm[m]][indexn[m]]=1;
    }
    double Z[K];
    for (int i=0;i<K;i++){
        Z[i]=s[i];
        for(int j=0;j<N;j++){
            if (S[i][j]==1){
                Z[i]=fabs(s[i]-S[i][j]*(j+1));
                break;
            }
        }
    }
    //check_int(K,N,D);
    double Y[K];
        for (int i=0;i<K;i++){
            int D_temp=0;
            for(int j=0;j<N;j++){
                D_temp=D_temp+D[i][j];
        }
        Y[i]=d[i]-D_temp;
    }


    int NN=4*l*N+l*l*N*M+l*N*M+2*K*N+2*K;
    double *x0=(double*)malloc( NN* sizeof(double*));
    for (int i=0; i<NN;i++){
        x0[i]=0;
    }
    
    for(int i=0;i<l;i++){
        for(int j=0;j<N;j++){
            x0[i+j*l]=(double) P[i][j];
            x0[l*N+i+j*l]=(double) L[i][j];
            x0[2*l*N+i+j*l]=(double) G[i][j];
            x0[3*l*N+l*l*N*M+l*N*M+2*K*N+i+j*l]=(double) H[i][j];
        }
    }
    
    for(int i=0;i<l;i++){
        for(int ii=0;ii<l;ii++){
            for(int im=0;im<M;im++){
                for(int j=0;j<N;j++){
                    if (F[i][ii][im][j]==1)
                        x0[3*l*N+l*l*M*j+l*l*im+ii*l+i]=(double) F[i][ii][im][j];
                }
            }
        }
    }
    
    for(int i=0;i<l;i++){
        for(int im=0;im<M;im++){
            for(int j=0;j<N;j++){
                if (E[i][im][j]==1)
                    x0[3*l*N+l*l*N*M+j*l*M+l*im+i]=E[i][im][j];
            }
        }
    }

    for(int k=0;k<K;k++){
        for (int j=0;j<N;j++){
            x0[3*l*N+l*l*N*M+l*M*N+k+j*K]=S[k][j];
            x0[3*l*N+l*l*N*M+l*M*N+K*N+k+j*K]=D[k][j];
        }
    }
    for(int k=0;k<K;k++){
        x0[4*l*N+l*l*N*M+l*M*N+2*K*N+k]=Z[k];
        x0[4*l*N+l*l*N*M+l*M*N+2*K*N+K+k]=Y[k];
    }
  //x=[P(:);L(:);G(:);F(:);E(:);S(:);D(:);H(:);Z;Y];
    return x0;
}
