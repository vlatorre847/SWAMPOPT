%N: Number of time intervals
%M: Number of operations
%l: Number of channels
%Ii: Set of the sets of the channels downstream every channel
%dt: Time interval duration in minutes
%tau: Delay times for every channel
%R0: Initial water stored in the channels
%H0: Initial opening of the gates
%V0: initial inlet flow in the channels
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

conf_file="synt_conf.txt";
req_file="synt_req.txt";

%conf_file="Real_World_conf.txt";
%req_file="Real_World_req.txt";

fp1 = fopen(conf_file, "r");
temp_l=fgetl(fp1);
N=sscanf(temp_l,"N=%f");
temp_l=fgetl(fp1);
M=sscanf(temp_l,"M=%f");
temp_l=fgetl(fp1);
l=sscanf(temp_l,"l=%f");
temp_l=fgetl(fp1);
[token,~]=strtok(temp_l,'Ii=');
Ii=sscanf(token,"%f",l);
temp_l=fgetl(fp1);
dt=sscanf(temp_l,"dt=%f");
temp_l=fgetl(fp1);
[token,remain]=strtok(temp_l,'tau=');
tau=sscanf(token,"%f",l);

delay=max(1,ceil(tau./dt));
delay_tot=sum(delay);

temp_l=fgetl(fp1);
[token,remain]=strtok(temp_l,'R_in=');
R0=sscanf(token,"%f",l);
temp_l=fgetl(fp1);
%For V0, the number of values depend on the delay on the channel. If for
%example the first canal has a delay of 2, and the second canal a delay of
%3 their values must be reported in the configuration file in the following
%fashion:
%P_in= 0 0, 0 0 0
%That is, the initial values for the channels are separated by a comma.
[token,remain]=strtok(temp_l,'V_in=');
V0=zeros(delay_tot,1);
for i=1:l
    for j=1:delay(i)
        V0(sum(delay(1:i-1))+j)=sscanf(token,"%f,");
    end
end
V0=V0*60*dt/1000;

temp=1;
V0_temp=cell(l,1);
for i=1:l
    V0_temp{i}=V0(temp:sum(delay(1:i)));
    temp=sum(delay(1:i))+1;
end
V0=V0_temp;


temp_l=fgetl(fp1);
[token,remain]=strtok(temp_l,'H_in=');
H0=sscanf(token,"%f",l);


temp_l=fgetl(fp1);
tp_size=sscanf(temp_l,"tp_size=%f");
temp_l=fgetl(fp1);
tpi_size=sscanf(temp_l,"tpi_size=%f");
temp_l=fgetl(fp1);
[token,remain]=strtok(temp_l,'tp=');
tp=sscanf(token,"%f",tp_size);
temp_l=fgetl(fp1);
[token,remain]=strtok(temp_l,'tpi=');
tpi=sscanf(token,"%f",tpi_size);
temp_l=fgetl(fp1);
[token,remain]=strtok(temp_l,'c=');
c=sscanf(token,"%f",l);
temp_l=fgetl(fp1);
[token,remain]=strtok(temp_l,'gamma=');
gamma=sscanf(token,"%f",l);
temp_l=fgetl(fp1);
[token,remain]=strtok(temp_l,'rho=');
rho=sscanf(token,"%f",l);
temp_l=fgetl(fp1);
[token,remain]=strtok(temp_l,'psi=');
psi1=sscanf(token,"%f",l*l);
temp_l=fgetl(fp1);
eps_min=sscanf(temp_l,"eps_min=10^%f");
eps_min=10^eps_min;
temp_l=fgetl(fp1);
[token,remain]=strtok(temp_l,'w=');
w=sscanf(token,"%f",4);


fclose(fp1);


c=c*60*dt/1000;
psi=zeros(l,l);
for i=1:l
    for j=1:l
        psi(i,j)=psi1(i+(j-1)*l)/dt;
    end
end
    


r=sum(c)*ones(N,1);
%delay=max(1,ceil(tau./dt));
eps_min=10^-8;


fp1 = fopen(req_file, "r");
temp_l=fgetl(fp1);
K=sscanf(temp_l,"K=%f");
temp_l=fgetl(fp1);
[token,~]=strtok(temp_l,'Ki=');
Ki=sscanf(token,"%f",K);
temp_l=fgetl(fp1);
[token,~]=strtok(temp_l,'q=');
q=sscanf(token,"%f",K);
temp_l=fgetl(fp1);
[token,~]=strtok(temp_l,'s=');
s=sscanf(token,"%f",K);
temp_l=fgetl(fp1);
[token,~]=strtok(temp_l,'d=');
d=sscanf(token,"%f",K);
temp_l=fgetl(fp1);
[token,~]=strtok(temp_l,'eps=');
eps=sscanf(token,"%f",K);
temp_l=fgetl(fp1);
[token,~]=strtok(temp_l,'alpha=');
alpha=sscanf(token,"%f",K);
temp_l=fgetl(fp1);
[token,~]=strtok(temp_l,'beta=');
beta=sscanf(token,"%f",K);

fclose(fp1);

q=q*60*dt/1000;
d=d*60/dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x0]=initial_solution(N,M,K,l,Ki,Ii,q,s,d,eps, gamma, psi, tp, tpi,c,rho,delay,H0,V0,R0);


Dt=max(max(s-1,N-s-eps.*d/dt),1);
Dv=(1-eps).*q.*d;


[f,A_beg,A_nnz,A_indx,A_val,A_sense,b,lb,ub,ctype,A_j,A_r]=get_MILP(N,M,K,l,tp_size,tpi_size,V0,R0,H0,tp,tpi,c,r,gamma,rho,psi,delay,eps_min,Ki,Ii,q,s,d,eps,alpha,beta,Dt,Dv,w);

 A=sparse(A_indx+1,A_j,A_val);
 [~,NN]=size(A);

%%%%%%intlinprog%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%intcon=zeros(l*M*N+2*K*N,1);
%k=1;
%for i=1:NN
%    if strcmp(ctype(i),'B')
%        intcon(k)=i;
%        k=k+1;
%    end
%end
%intcon(intcon==0)=[];
%options=optimoptions(@intlinprog);
%options.MaxTime=3600;
%[x,fval,exitflag,output] = intlinprog(f,intcon,A(1:A_r,:),b(1:A_r),A(A_r+1:end,:),b(A_r+1:end),lb,ub,x0,options);
%%%%%%intlinprog%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%CPLEX%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   options = cplexoptimset('cplex');
   options.timelimit=3600;
   options.display='on';
   [x,fval,flag,output] = cplexmilp(f,A(1:A_r,:),b(1:A_r),A(A_r+1:end,:),b(A_r+1:end),[],[],[],lb,ub,ctype,x0,options);
%%%%%%CPLEX%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V=reshape(x(1:N*l),l,N);
R=reshape(x(N*l+1:2*N*l),l,N);
G=reshape(x(2*N*l+1:3*N*l),l,N);
F=reshape(x(3*N*l+1:3*N*l+l*l*N*M),l,l,M,N);
E=reshape(x(3*N*l+l*l*N*M+1:3*N*l+l*l*N*M+l*N*M,1),l,M,N);
S=reshape(x(3*N*l+l*l*N*M+l*N*M+1:3*N*l+l*l*N*M+l*N*M+K*N),K,N);
D=reshape(x(3*N*l+l*l*N*M+l*N*M+K*N+1:3*N*l+l*l*N*M+l*N*M+K*N+K*N),K,N);
H=reshape(x(3*N*l+l*l*N*M+l*N*M+K*N+K*N+1:3*N*l+l*l*N*M+l*N*M+K*N+K*N+l*N,1),l,N);
