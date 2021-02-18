function [x]=initial_solution(N,M,K,l,Ki,Ii,q,s,d,eps, gamma, psi, tp, tpi,c,rho,delay,H0,V0,R0)
%N: Number of time intervals
%M: Number of operations
%K: Numer of irrigations
%l: Number of channels
%Ki: Set of the sets of off-takes on the channels
%Ii: Set of the sets of the channels downstream every channel
%q: Quantity of water required by the off-take for time interval
%s: Desidered starting time interval for the irrigation
%d: Desidered duration for the irrigation 
%eps: minimum ration of water that must be delivered for irrigation
%gamma: seepage ration for channel
%psi: travelling time from one gate to another
%tp: Time intervals where the gate-keeper cannot operate
%c: Maximum inlet volume capacity for every gate
%rho: Minimum ratio of water that must flow through an open gate
%delay: time intervals the water takes to go through a channel
%R0: Initial water stored in the channels
%H0: Initial opening of the gates
%V0: initial inlet flow in the channels

P=zeros(l,N);
L=P;
G=P;
H=P;
F=zeros(l,l,M,N);
E=zeros(l,M,N);
S=zeros(K,N);
D=S;
dappo=P; % How many irrigation are active for every channel in each time interval
dmin=ceil(d.*eps); %Minimum duration for irrigation
prec=Ii;   %upstream channel for every channel

Ii_temp=cell(l,1);
for i=2:l
    temp=Ii(i);
    Ii_temp{temp}=[Ii_temp{temp};i];
end
Ii=Ii_temp;



Ki_temp=cell(l,1);
for i=1:K
    temp=Ki(i);
    Ki_temp{temp}=[Ki_temp{temp};i];
end
Ki=Ki_temp;


%%%%%%%%Initial Conditions on the channels channel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%P(:,1)=V0;


%%%%%%%%P,L,S,D in the first (source) channel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startj=zeros(l,1); % Time interval in which the channel has been opened
startj(1)=ceil(psi(1,1));
G(1,startj(1))=1;

P(1,startj(1):N)=c(1);

for ik=1:length(Ki{1})
    j=max(delay(1)+startj(1),s(Ki{1}(ik)));
    if j+dmin(Ki{1}(ik))>N  %If you cannot do the irrigation in the considered time interval break
        continue;
    end

    while j<=N
        if P(1,j-1)+L(1,j-1)-sum(q(Ki{1}).*D(Ki{1},j-1+delay(1)))<q(ik)
            j=j+1;
            continue
        end
        S(Ki{1}(ik),j)=1;
        D(Ki{1}(ik),j:j+dmin(Ki{1}(ik))-1)=1;
        break
    end
end
stopj=ones(l,1)*N+1;
j=2;
temppsi=0;
%For cycle for opening the channels, so that the time necessary to do the
%operations does not exceede the duration of the time interval itself
for i=2:l
    while j<N
        if temppsi+psi(i-1,i)>1
            j=j+1;
            temppsi=temppsi-1;%temppsi-1;
            continue
        end
        startj(i)=j;
        %G(i,j)=1;
        if psi(i-1,i)>1
            temppsi=0;
        else
            temppsi=temppsi+psi(i-1,i);
        end
        break;
    end
end
H(1,1:N)=P(1,1:N)/c(1);
for i=2:l
    for j=1:N
        P(i,j)=V0{i}(delay(i));
        H(i,j)=H0(i);
    end
end


%%%%%%%%S, D, G%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[L,P,G,H]= add_irrigation(L,P,G,H,D,stopj,startj,q,c,Ii,Ki,gamma,psi,l,N,tp,rho,prec,delay,V0,R0);%Update the flows after adding an irrigation        

for ii=2:l        
    [~,indexk]=sort(dmin(Ki{ii}),'descend');
    for ik=1:length(Ki{ii})
        flag=0;
        t=indexk(ik);
        j=max(startj(ii)+1,s(Ki{ii}(t)));
        while j<=N
            if j+dmin(Ki{ii}(t))>N&&flag==0
                flag=1;
                j=startj(ii)+1;
            end
            if j+dmin(Ki{ii}(t))>N &&flag==1 %If you cannot do the irrigation in the considered time interval break
                break;
            end
            if ismember(j,tpi)
                  j=j+1;
                 continue
            end

            D1=D;
            D1(Ki{ii}(t),j:j+dmin(Ki{ii}(t))-1)=1;
            [L1,P1,G1,H1]= add_irrigation(L,P,G,H,D1,stopj,startj,q,c,Ii,Ki,gamma,psi,l,N,tp,rho,prec,delay,V0,R0);% add an irrigation at time j  
            
            if (min(min(L1))<-10^-10||min(min(P1))<-10^-10) % Is this irrigation possible?
                  j=j+1;
                 continue
            end
               
            
            S(Ki{ii}(t),j)=1;
            D=D1;
            L=L1;
            P=P1;
            G=G1;
            H=H1;
            % Now the irrigation has been added
            break
        end

        for i=1:l
            for j=1:N
                dappo(i,j)=sum(D(Ki{i},j));
            end
        end       
    end
    temp=find(dappo(ii,:)>0,1,'last');
    if ~isempty(temp)
        while ismember(temp,tp)
            temp=temp+1;
        end
        stopj(ii)=temp; % Time interval when the last irrigation has been performed in the channel, and it can be closed 
    end
end

%%%%%%%%E%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im=nnz(G);        
m=M;
k=0;
indexi=zeros(im,1);
indexm=zeros(im,1);
indexn=zeros(im,1);
for n=N:-1:1
    for i=l:-1:1
        if G(i,n)==1
            E(i,m,n)=1;
            indexi(im-k)=i;
            indexm(im-k)=m;
            indexn(im-k)=n;
            k=k+1;
            m=m-1;
        end
    end
end
%%%%%%%%F%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F(indexi(1),indexi(1),indexm(1),indexn(1))=1;
for m=2:im
    F(indexi(m),indexi(m-1),indexm(m),indexn(m))=1;
end

Z=abs(s-S*(1:N)');
Y=d-sum(D,2);
x=[P(:);L(:);G(:);F(:);E(:);S(:);D(:);H(:);Z;Y];
index_rho=sum(D,2)==0;
eps(index_rho)=0;


end



function [R,V,G,H]= add_irrigation(R,V,G,H,D,stopj,startj,q,c,Ii,Ki,gamma,psi,l,N,tp,rho,prec,delay,V0,R0)
%This function Adjusts the flows on the network according to the water
%balance constrain after the new irrigation has been temporarily added to
%the scheduling. In the main program, if just one of the V or R are
%negative, this new scheduling is discarted.
V(1,startj(1):N)=c(1);
CR=chan_requ(D,N,l,gamma,q,prec,Ki,tp,delay);
CR=[CR,zeros(l,1)];
for i=1:l
    for j=1:delay(i)
        if j==1
            L_temp=R0(i);
        else
            L_temp=R(i,j-1);
        end
        for k=1:length(Ii{i})
            ii=Ii{i}(k);
            psi_temp=gk_workload(j,ii,psi,G,N);
            tempf=((1-gamma(i))^delay(i))*V0{i}(j)+(1-gamma(i))*L_temp - sum(q(Ki{i}).*D(Ki{i},j))-sum(CR(Ii{i}(Ii{i}~=ii),j));
            if ~isempty(Ii{ii})
                V(ii,j)=min([tempf,c(ii)]);
            else
                tempd=sum(q(Ki{ii}).*D(Ki{ii},min(j+delay(ii),N)))/((1-gamma(ii))^delay(ii));
                V(ii,j)=min([tempf,c(ii),tempd]);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Set the values for H
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [V,G,H,CR]=adjust_flow(V,G,H,CR,ii,j,rho,psi,N,tp,c,psi_temp,V0,delay);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Set the values for R
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        R(i,j)=((1-gamma(i))^delay(i))* V0{i}(j)        +(1-gamma(i))*L_temp  -sum(q(Ki{i}).*D(Ki{i},j))-sum(V(Ii{i},j));

    end
    
    
    for j=1+delay(i):N
        for k=1:length(Ii{i})
            ii=Ii{i}(k);
            psi_temp=gk_workload(j,ii,psi,G,N);
            tempf=((1-gamma(i))^delay(i))*(V(i,j-delay(i)))+(1-gamma(i))*R(i,j-1) - sum(q(Ki{i}).*D(Ki{i},j))-sum(CR(Ii{i}(Ii{i}~=ii),j));

            if ~isempty(Ii{ii})
                V(ii,j)=min([tempf,c(ii)]);

            else
                tempd=sum(q(Ki{ii}).*D(Ki{ii},min(j+delay(ii),N)))/((1-gamma(ii))^delay(ii));
                V(ii,j)=min([tempf,c(ii),tempd]);
                if j>=stopj(ii)||j<startj(ii)
                    V(ii,j)=0;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Set the values for H
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [V,G,H,CR]=adjust_flow(V,G,H,CR,ii,j,rho,psi,N,tp,c,psi_temp,V0,delay);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Set the values for R
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        R(i,j)=((1-gamma(i))^delay(i))*(V(i,j-delay(i)))+(1-gamma(i))*R(i,j-1)-sum(q(Ki{i}).*D(Ki{i},j))-sum(V(Ii{i},j));
    end
    for ii=2:l
        for jj=2:N
            if H(ii,jj)==H(ii,jj-1)&&G(ii,jj)==1
                G(ii,jj)=0;
            end
        end
    end
end
end


function CR= chan_requ(D,N,l,gamma,q,prec,Ki,tp,delay)
%Computes the water required for every channel in every time interval
CR=zeros(l,N);
        for i=1:l
            for j=1+delay(i):N
                if ismember(j-delay(i),tp)
                    CR(i,j-delay(i))=CR(i,j-delay(i)-1);
                else
                    CR(i,j-delay(i))=sum(q(Ki{i}).*D(Ki{i},j))/(1-gamma(i))^(delay(i));
                end
            end
        end
        for i=1:l
            temp=delay(i);
            prec_c=prec(i);
            while prec_c~=0
                temp=temp+delay(prec_c);
                prec_c=prec(prec_c);
            end
            for j=temp:N
                if CR(i,j)==0
                    continue
                end
                prec_c=prec(i);
                jj=0;
                appo=CR(i,j);
                while prec_c~=0
                    jj=jj+delay(prec_c);
                    appo=appo/(1-gamma(prec_c))^jj;
                    CR(prec_c,j-jj)=CR(prec_c,j-jj)+appo;
                    prec_c=prec(prec_c);
                end
            end
        end
end


function psi_temp=gk_workload(j,ii,psi,G,N)
%Checks if the gatekeeper can perform the operation on gate ii on time j
last_gate=[];
if j>1
    jj=1;
    while isempty(last_gate)
        last_gate=find(G(:,j-jj),1,'last');
        jj=jj+1;
    end
    index1=find(G(:,j)==1);
    if isempty(index1) %If there is no operation in time interval j
        if psi(last_gate,ii)<=1
            psi_temp=psi(last_gate,ii);
        else
            psi_temp=-jj+1+psi(last_gate,ii)+1;
        end
    else
        G_appo=G(:,j);
        G_appo(ii)=1;
        index1=find(G_appo==1);
        if psi(last_gate,index1(1))<=1
            psi_temp=psi(last_gate,index1(1));
        else
            %psi_temp=j-jj+psi(last_gate,ii);
            psi_temp=0;
        end
        for iii=2:length(index1)
            psi_temp=psi_temp+psi(index1(iii-1),index1(iii));
        end
    end
else
    index1=find(G(:,j)==1);
    if isempty(index1) %If there is no operation in time interval j
        psi_temp=0;
    else
        G_appo=G(:,j);
        G_appo(ii)=1;
        index1=find(G_appo==1);
        psi_temp=psi(index1(1),index1(1));
        for iii=2:length(index1)
            psi_temp=psi_temp+psi(index1(iii-1),index1(iii));
        end
    end
end

jn=j;
while jn~=N %Check if the new operation does not interfere with the future operations
    if sum(G(:,jn+1))==0
        jn=jn+1;
        continue
    end
    G_appo1=G(:,jn+1);
    index2=find(G_appo1==1);
    if psi(ii,index2(1))<=1
        psi_temp1=psi(ii,index2(1));
        for iii=2:length(index2)
            psi_temp1=psi_temp1+psi(index2(iii-1),index2(iii));
        end
    else
        psi_temp1=j-jn+psi(ii,index2(1));
    end
    psi_temp=max(psi_temp,psi_temp1);
    break
end
end

function [V1,G1,H1,CR1]=adjust_flow(V,G,H,CR,ii,j,rho,psi,N,tp,c,psi_temp,V0,delay)
%Checks if it is really possible to adjust the flow to the suggested value
%in V for channel ii in time interval j, by checking if the gatekeeper can perform the necessary operations.
%If it is not possible, leaves the value of V equal to that of the previous time interval.
%It modifies G and H accordingly.
if j<=delay(ii)
    P_temp=V0{ii}(j);
else
    P_temp=V(ii,j-1);
end

if V(ii,j)<rho(ii)*H(ii,j)*c(ii)
    if (ismember(j,tp)||psi_temp>1)
        V(ii,j)=P_temp;
        CR(ii,j)=max(V(ii,j),CR(ii,j));
        V1=V;
        G1=G;
        H1=H;
        CR1=CR;
        return;
    else
        G(ii,j)=1;
        jj=j+1;
        while jj<N+1
            psi_temp1=gk_workload(jj,ii,psi,G,N);
            if psi_temp1<=1 || G(ii,jj)==1
                if (ismember(jj,tp))
                    jj=jj+1;
                    continue
                end
                G(ii,jj)=1;
                break
            end
            jj=jj+1;
        end
        H(ii,j:jj-1)=min(V(ii,j:jj-1))/(c(ii)*rho(ii));
    end
else
    if V(ii,j)>H(ii,j)*c(ii)
        if (ismember(j,tp)||psi_temp>1)
            V(ii,j)=P_temp;
            CR(ii,j)=max(V(ii,j),CR(ii,j));
            V1=V;
            G1=G;
            H1=H;
            CR1=CR;
            return;
        else
            G(ii,j)=1;
            jj=j+1;
            while jj<N+1
                psi_temp1=gk_workload(jj,ii,psi,G,N);
                if psi_temp1<=1 ||G(ii,jj)==1
                    if (ismember(jj,tp))
                        jj=jj+1;
                        continue
                    end
                    G(ii,jj)=1;
                    break
                end
                jj=jj+1;
            end
            H(ii,j:jj-1)=max(V(ii,j:jj-1))/c(ii);
        end
    end
end
V1=V;
G1=G;
H1=H;
CR1=CR;
return;

end
