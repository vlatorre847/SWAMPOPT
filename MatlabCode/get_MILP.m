function [f,A_beg,A_nnz,A_indx,A_val,A_sense,b,lb,ub,ctype,A_j, A_r]=get_MILP(N,M,K,l,tp_size,tpi_size,V0,R0,H0,tp,tpi,c,r,gamma,rho,psi,delay,eps_min,Ki,Ii,q,s,d,eps,alfa,beta,Dt,Dv,j_funct)
% global A_11	A_12	A_21	A_22	A_3	A_4	A_5	A_6	A_7	A_8	A_9	A_101	A_102 
% global Aeq_1 Aeq_2 Aeq_3 Aeq_4 Aeq_6 Aeq_7 Aeq_8 A_2eqini 
%N: Number of time intervals
%M: Number of operations
%l: Number of channels
%Ii: Set of the sets of the channels downstream every channel
%dt: Time interval duration in minutes
%tau: Delay times for every channel
%R0: Initial water stored in the channels
%H0: Initial opening of the gates
%tp_size: number of rest interval for the gatekeeper
%tpi_size: number of intervals it is not possible to start irrigations
%tp: Time intervals where the gate-keeper cannot operate
%tpi: Time intervals the irrigations cannot start
%c: Maximum inlet volume capacity for every gate
%rho: Minimum ratio of water that must flow through an open gate
%psi: travelling time from one gate to another
%w: weights of the objective function

%K: Numer of irrigations
%Ki: Set of the sets of off-takes on the channels
%q: Quantity of water required by the off-take per time interval
%s: Desidered starting time interval for the irrigation
%d: Desidered duration for the irrigation 
%eps: minimum ratio of water that must be delivered for irrigation
%alfa: time priority for the k-th irrigation
%beta: volume priority for the k-th irrigation
%Dt: Maximum delay for an irrigation


%delay: time intervals the water takes to go through a channel
%Dv: Minimum volume that can be delivered
%Constraints creation, latex notation used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Aeq1
%(1-\gamma_i)^{\tau_i}V_i^{n-\tau_i} +(1-\gamma_i)R_i^{n-1}- \sum_{k\in K_i} q_kD_k^n-\sum_{j\in I_i}V_j^n-R_i^n=0
appo=zeros(l);
appo1=zeros(l,K);
temp=zeros(l,(max(delay)+1)*l);
templ=zeros(l,(max(delay)+1)*l);
temp1=zeros(l*N,(N+max(delay))*l);
templ1=zeros(l*N,(N+max(delay))*l);
tempd=zeros(l,K*(max(delay)+1));
tempd1=zeros(l*N,K*(N+max(delay)));
for i=1:l
    appo(i,Ii==i)=-1;
    for kk=1:K
        if Ki(kk)==i
            appo1(i,kk)=-q(kk);
        end
    end
   temp(i,(delay(i))*l+1 :(delay(i)+1)*l)=appo(i,:);
   templ(i,delay(i)*l+i)=1;
   tempd(i,(delay(i))*K+1:(delay(i)+1)*K)=appo1(i,:);
end
for j=1:N
   temp1((j-1)*l+1:j*l,(j-1)*l+1:(j-1)*l+max((delay)+1)*l)= temp;
   templ1((j-1)*l+1:j*l,(j-1)*l+1:(j-1)*l+max((delay)+1)*l)= templ;
   tempd1((j-1)*l+1:j*l,(j-1)*K+1:(j-1)*K+max((delay)+1)*K)=tempd;
end

temp1(:,N*l+1:end)=[];
templ1(:,N*l+1:end)=[];
tempd1(:,K*N+1:end)=[];

ini_P=zeros(l*(max(delay)-1),l*(max(delay)-1));
ini_D=zeros(l*(max(delay)-1),K*(max(delay)-1));
ini_L=zeros(l*(max(delay)-1),l*(max(delay)));
appo_L=eye(l);
delay_vec=zeros(l,max(delay)-1);
delay_mat=zeros(max(delay)-1,max(delay)-1,l);

for i=1:l
    delay_vec(i,1:delay(i)-1)=1;
    delay_mat(:,:,i)=diag(delay_vec(i,:));
    ini_P((i-1)*(max(delay)-1)+1:i*(max(delay)-1),:)= kron(delay_mat(:,:,i),appo(i,:));
    ini_D((i-1)*(max(delay)-1)+1:i*(max(delay)-1),:)= kron(delay_mat(:,:,i),appo1(i,:));
    ini_L((i-1)*(max(delay)-1)+1:i*(max(delay)-1),:)= [kron(delay_mat(:,:,i)*(1-gamma(i)),appo_L(i,:)),zeros(max(delay)-1,l)]-[zeros(max(delay)-1,l),kron(delay_mat(:,:,i),appo_L(i,:))];
end

index_ini=delay_vec';
index_ini=find(index_ini==0);

ini_P(index_ini,:)=[];
ini_D(index_ini,:)=[];
ini_L(index_ini,:)=[];
[temp_size,~]=size(ini_P);

temp1=temp1+kron(spdiags(ones(N,1),0,N,N),diag((1-gamma)).^delay);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Aeqini and Aeqini1

%\sum_{k\in K_i} q_kD_k^1+\sum_{j\in I_i}V_j^1+(1-\gamma_i)(\bar{R}_i^0+\bar{V}_i^0)= R^1_i
%\sum_{k\in K_i} q_kD_k^n+\sum_{j\in I_i}V_j^n+(1-\gamma_i)(R_i^{n-1}+\bar{V}_i^{n-1})= R^{n}_i

Aeqini=[appo,sparse(l,l*(N-1)),-eye(l),sparse(l,l*(N-1)),sparse(l,N*l+l*l*N*M+l*N*M),sparse(l,N*K),appo1,sparse(l,(N-1)*K),sparse(l,N*l)];
Aeqini1=[sparse(temp_size,l),ini_P,sparse(temp_size,l*(N-max(delay))),ini_L,sparse(temp_size,(N-(max(delay)))*l),sparse(temp_size,N*l+l*l*N*M+l*N*M+K*N),sparse(temp_size,K),ini_D,sparse(temp_size,K*(N-max(delay))),sparse(temp_size,N*l)];
templ2=templ1(:,l+1:l*N);
appo_rep=repmat((1-gamma),N,1);
templ2=[templ2,zeros(N*l,l)];
for i=1:l*N
templ2(:,i)=templ2(:,i).*appo_rep;
end

Aeq1=[temp1,templ2-templ1,sparse(N*l,N*l),sparse(N*l,l*l*N*M),sparse(N*l,l*N*M),sparse(N*l,K*N),tempd1,sparse(N*l,N*l)];
delay_temp=N-delay;
index_canc=[];
for i=1:l
   index_canc=[index_canc;((delay_temp(i):N-1)*l+i)'];
end
index_canc=sort(index_canc);

 for i=sum(delay):-1:1
     aeq1_temp1=Aeq1(1:index_canc(i)-1,:);
     aeq1_temp2=Aeq1(index_canc(i)+1:end,:);
     Aeq1=[aeq1_temp1;aeq1_temp2];
 end

Aeq1=[Aeqini;Aeqini1;Aeq1];
%Initial flows must be with a minus in front of them

beqini=zeros(l,1);
for i=1:l
    beqini(i)=-(1-gamma(i)).*(R0(i))-(1-gamma(i)).^(delay(i))*V0{i}(1);
end

beqini1=zeros(temp_size,1);
for i=1:l
    for j=2:delay(i)
        beqini1(j-1)=-((1-gamma(i))^(delay(i))).*V0{i}(1);
    end
end
beq1=zeros(l*N-sum(delay),1);
beq1=[beqini;beqini1;beq1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A11
%H_i^n-H_i^{n-1}\le \Omega_i G_i^n
appo=-spdiags([ones(N-1,1);0],0,N,N)+spdiags(ones(N,1),1,N,N);
A11=[sparse(N*l,N*l),sparse(N*l,N*l),-diag(repmat(ones(l,1),N-1,1),l),sparse(N*l,l*l*N*M+l*N*M+2*K*N),kron(appo,spdiags(ones(l,1),0,l,l))];
%A12
%H_i^n-H_i^{n-1}\ge-\Omega_i G_i^n
A12=[sparse(N*l,N*l),sparse(N*l,N*l),-diag(repmat(ones(l,1),N-1,1),l),sparse(N*l,l*l*N*M+l*N*M+2*K*N),-kron(appo,spdiags(ones(l,1),0,l,l))];
A11=A11(1:(N-1)*l,:);
A12=A12(1:(N-1)*l,:);

b11=ones(l*(N-1),1)*eps_min;
b12=-ones(l*(N-1),1)*eps_min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Aini1
%H_i^1-\Omega_i G_i^1&\le H_i^{0} 
%-H_i^1-\Omega_i G_i^1&\le H_i^{0}
%initial gate opening conditions: The gate opening ratio cannot change at the beginning if the gate is not operated 
Aini11=[sparse(l,l*N),sparse(l,l*N),-eye(l),sparse(l,(N-1)*l),sparse(l,l*l*N*M+l*N*M+2*K*N),eye(l),sparse(l,(N-1)*l)];
Aini12=[sparse(l,l*N),sparse(l,l*N),-eye(l),sparse(l,(N-1)*l),sparse(l,l*l*N*M+l*N*M+2*K*N),-eye(l),sparse(l,(N-1)*l)];
bini1=H0;
bini2=-bini1;

A11=[Aini11;A11];
A12=[Aini12;A12];
b11=[bini1;b11];
b12=[bini2;b12];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A2
%V_i^n\le  min\{c_i,r_n \} 
%-V_i^n+H_i^n \rho_i \min \{c_i,r_n\}\le 0
A21=[spdiags(ones(N*l,1),0,N*l,N*l),sparse(N*l,2*N*l+l*l*N*M+l*N*M+2*K*N),-kron(spdiags(ones(N,1),0,N,N),diag(c))];
b21=zeros(N*l,1);

A22=[-spdiags(ones(N*l,1),0,N*l,N*l),sparse(N*l,2*N*l+l*l*N*M+l*N*M+2*K*N),kron(spdiags(ones(N,1),0,N,N),diag(c.*rho))];
b22=zeros(N*l,1);



A2eqini=[kron(eye(N),[1,zeros(1,l-1)]),sparse(N,2*l*N+l*l*N*M+l*N*M+2*K*N),-kron(eye(N),[c(1),zeros(1,l-1)])];
b2eqini=zeros(N,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Aeq2
%\sum_{i=1}^l\ G_i^n= 0, for n \in t^p

Aeq2=zeros(tp_size,3*l*N+l^2*N*M+l*N*M+2*K*N+l*N);% for G

for i=1:tp_size
   Aeq2(i,2*l*N+(tp(i)-1)*l+1:2*l*N+(tp(i))*l)=ones(l,1); 
end
beq2=zeros(tp_size,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A3
% \sum_{i=1}^{l}\sum_{n=1}^N E_i^{n,m}\le 1  
A3=[sparse(M,3*l*N+N*M*l*l),repmat(kron(spdiags(ones(M,1),0,M,M),ones(1,l)),1,N),sparse(M,2*K*N),sparse(M,l*N)];
b3=ones(M,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Aeq3
%G_i^n- \sum_{m=1}^M E_i^{n,m}&=0
Aeq3=[sparse(l*N,2*l*N),eye(l*N),sparse(l*N,l^2*N*M),-kron(spdiags(ones(N,1),0,N,N),repmat(spdiags(ones(l,1),0,l,l),1,M)),sparse(l*N,2*K*N),sparse(l*N,l*N)];
beq3=zeros(l*N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A4
%E_i^{m-1,n}-\sum_{j=1}^l\sum_{p=n}^N F_{ij}^{m,p}\le0
appo=diag(ones((M-1)*l,1),l);
appo((M-1)*l+1:M*l,:)=[];
appo=sparse(appo);
appo=-kron(kron(sparse(triu(ones(N))),appo),ones(1,l));       
appo1=eye(l*M);
appo1((M-1)*l+1:M*l,:)=[];
appo1=sparse(appo1);
A4=[sparse(l*N*(M-1),3*l*N),appo,kron(spdiags(ones(N,1),0,N,N),appo1),sparse(l*N*(M-1),2*K*N),sparse(l*N*(M-1),l*N)];
b4=zeros(l*N*(M-1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Aeq4
%E_j^{m,n}= \sum_{i=1}^l F_{ij}^{m,n}
Aeq4=[sparse(l*N*M,3*l*N),-kron(spdiags(ones(M*N,1),0,M*N,M*N),repmat(spdiags(ones(l,1),0,l,l),1,l)),spdiags(ones(l*M*N,1),0,l*M*N,l*M*N),sparse(l*N*M,2*K*N),sparse(l*N*M,l*N)];
beq4=zeros(l*N*M,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A5
%\sum_{j=1}^l\sum_{n=1}^N nE_j^{n,m}-\sum_{j=1}^l\sum_{n=1}^N nE_j^{n,m-1}\ge  \sum_{i=1}^l\sum_{j=1}^l\sum_{n=1}^N \psi_{ij}F_{ij}^{n,m}
psi_temp=psi(:);
psi_temp(psi_temp<=1)=0;
A5=[sparse(M,3*l*N),kron(repmat(spdiags(ones(M,1),1,M,M),1,N),psi_temp'),1*kron(kron(1:N,spdiags(ones(M,1),0,M,M)-spdiags(ones(M,1),1,M,M)),ones(1,l)),sparse(M,2*K*N),sparse(M,l*N)];
A5(M,:)=[];
b5=zeros(M-1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A6
% \sum_{i=1}^l\sum_{j=1}^l\sum_{m=1}^M \psi_{ij}F_{ij}^{n,m}\le 1
psi_temp=psi(:);
psi_temp(psi_temp>1)=0;

A6=[sparse(N,3*l*N),kron(spdiags(ones(N,1),0,N,N),repmat(psi_temp',1,M)),sparse(N,N*M*l),sparse(N,2*K*N),sparse(N,l*N)];
b6=ones(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A7
% \sum_{n=1}^N S_k^n\le1
A7=[sparse(K,3*l*N+N*M*l+N*M*l*l),repmat(spdiags(ones(K,1),0,K,K),1,N),sparse(K,K*N),sparse(K,l*N)];
b7=ones(K,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Aeq5
% S_k^1= D_k^1
appo=[eye(K),sparse(K,K*(N-1))];
Aeq5=[sparse(K,3*l*N+N*M*l+N*M*l*l),appo,-appo,sparse(K,l*N)];
beq5=zeros(K,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A8
% S_k^n=0\rightarrow D_k^{n-1}-D_k^n\ge0
appo=-diag(ones(K*(N-1),1),K);
appo(K*(N-1)+1:K*N,:)=[];
appo1=eye(K*N)-diag(ones(K*(N-1),1),-K);
appo1(1:K,:)=[];

A8=[sparse(K*(N-1),3*l*N+N*M*l+N*M*l*l),appo,appo1,sparse(K*(N-1),l*N)];
b8=zeros(K*(N-1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A9
% \sum_{n=1}^N nS_k^n+\sum_{n=1}^N D_k^n&\le N
A9=[sparse(K,3*l*N+N*M*l+N*M*l*l),kron(1:N,spdiags(ones(K,1),0,K,K)),kron(ones(1,N),spdiags(ones(K,1),0,K,K)),sparse(K,l*N)];
b9=ones(K,1)*N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A10
%\epsilon_k q_k d_k\sum_{n=1}^N D_k^n\le q_k\sum_{n=1}^N  D_k^n\le q_k d_k 
 A101=[sparse(K,3*l*N+N*M*l+N*M*l*l+K*N),kron(ones(1,N),spdiags(q,0,K,K)),sparse(K,l*N)];
 b101=q.*d;

A102=-[sparse(K,3*l*N+N*M*l+N*M*l*l),repmat(spdiags(-eps.*d.*q,0,K,K),1,N),kron(ones(1,N),spdiags(q,0,K,K)),sparse(K,l*N)];
b102=zeros(K,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Aeq6
%\sum_{k=1}^K\ S_k^n= 0, for n \in t^{pi}

Aeq6=zeros(tpi_size,3*l*N+l^2*N*M+l*N*M+2*K*N+l*N);% for S

for i=1:tpi_size
    Aeq6(i,3*l*N+l^2*N*M+l*N*M+(tpi(i)-1)*K+1:3*l*N+l^2*N*M+l*N*M+(tpi(i))*K)=ones(K,1);
end
beq6=zeros(tpi_size,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Az1  -z_k-\sum_{n=1}^N n S_k^n\le -s_k  \forall k=1,\dots,K\\
Az1=[sparse(K,3*l*N+N*M*l*l+l*M*N),-kron(1:N,eye(K)),sparse(K,K*N+l*N),-eye(K)];
%Az2  -z_k+\sum_{n=1}^N n S_k^n\le s_k  \forall k=1,\dots,K\\
Az2=[sparse(K,3*l*N+N*M*l*l+l*M*N),kron(1:N,eye(K)),sparse(K,K*N+l*N),-eye(K)];
bz1=-s;
bz2=s;

%Ay  y_k+\sum_{n=1}^N D_k^n=d_k \forall k=1,\dots,K\\
Ay=[sparse(K,3*l*N+N*M*l*l+l*M*N+K*N),kron(ones(1,N),eye(K)),sparse(K,K+l*N),eye(K)];
by=d;



PSI=N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%A8=[sparse(N*l,N*l),eye(N*l),sparse(N*l,2*N*l+l*l*N*M+l*N*M+2*K*N)];
%b8=repmat(ubL,[N,1]);




lb=zeros(3*l*N+l^2*N*M+l*N*M+2*K*N+K+K+l*N,1);
ub=[max(r)*ones(2*l*N,1);ones(l*N+l^2*N*M+l*N*M+2*K*N+l*N,1);N*ones(K,1);d];


%%%%
A=[A11;A12;A21;A22;A3;A4;A5;A6;A7;A8;A9;A101;A102];
b=[b11;b12;b21;b22;b3;b4;b5;b6;b7;b8;b9;b101;b102];
Aeq=[Aeq1;Aeq2;Aeq3;Aeq4;Aeq5;Aeq6;A2eqini];
beq=[beq1;beq2;beq3;beq4;beq5;beq6;b2eqini];
%%%%


%NN=4*l*N+l^2*N*M+l*N*M+2*K*N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[s1,~]=size(A);
A=[A            ,sparse(s1,K);
   [Az1;Az2]];
[s1,~]=size(Aeq);
Aeq=[Aeq,sparse(s1,K)];
b=[b;bz1;bz2];
%NN=NN+K;

[s1,~]=size(A);

A=[A,sparse(s1,K)];
[s1,~]=size(Aeq);
Aeq=[Aeq,sparse(s1,K);
     Ay];
beq=[beq;by];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j=j_funct;
appo=j(1)+j(2)+j(3)+j(4);
j(1)=j(1)/appo;
j(2)=j(2)/appo;
j(3)=j(3)/appo;
j(4)=j(4)/appo;
psi_temp=psi(:);

f=[zeros(l*N,1); j(3)*1/(sum(r))*ones(l*N,1); zeros(l*N,1);j(4)*kron(ones(N*M,1),psi_temp)/(PSI); zeros(l*N*M,1) ; zeros(N*K,1) ; zeros(N*K,1);zeros(N*l,1);j(1)*alfa.*ones(K,1)/(sum(Dt));j(2)*beta.*ones(K,1)/(sum(Dv))];


ctypePL=repelem('C',N*l);
ctypeG=repelem('C',N*l);
ctypeF=repelem('C',l*l*M*N);
ctypeE=repelem('B',l*M*N);
ctypeS=repelem('B',K*N);
ctypeD=repelem('B',K*N);
ctypeZ=repelem('C',K);




for i=1:tp_size
    ctypeE(l*(tp(i)-1)*M+1:l*tp(i)*M)=repelem('C',l*M);
end



ctype=sprintf('%s%s%s%s%s%s%s%s%s%s',ctypePL,ctypePL,ctypeG,ctypeF,ctypeE,ctypeS,ctypeD,ctypePL,ctypeZ,ctypeZ);

[A_r,n_c]=size(A);
[Aeq_r,~]=size(Aeq);


sense_A=repelem('L',A_r);
sense_Aeq=repelem('E',Aeq_r);
A_sense=sprintf('%s%s',sense_A,sense_Aeq);

A=[A;Aeq];
b=[b;beq];
A_nnz=full(sum(A~=0)');
[A_indx,A_j]=find(A);
A_indx=A_indx-1;
A_val=nonzeros(A);

A_beg=zeros(n_c,1);
temp_sum=0;
for i=1:n_c
    A_beg(i)=temp_sum;
    temp_sum=temp_sum+A_nnz(i);
end


end
