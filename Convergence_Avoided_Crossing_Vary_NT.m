clear all
close all

%==========================================================================
%   Initalise creation and anhllation Operators 
%==========================================================================

N=128
num_efn=2;
hbar_eff=0.2;
dN=1
F=0.3
Qmax=1.5;
Pmax=1.5;
alpha=1;
gamma=0.001
beta=0;
set_efn='G'; % Invarient Subspace: Gain ('G') Loss ('L')
% Rescaling
alpha=alpha*hbar_eff^2;
F=F/hbar_eff;
Qmax=Qmax/hbar_eff;
Pmax=Pmax/hbar_eff;
[a,Q,P]=init_number_basis(N,1,1,1); % Get the operators in the number basis
H0=0.5*(P^2)-0.5*beta*Q^2+0.25*alpha*(Q^4)-0.5*1i*gamma*P^2;
NTmin=400;
NTmax=600;
dNT=3;
NTR=linspace(NTmin,NTmax,dNT)
% return
tic
for itt=1:length(NTR)
NT=NTR(itt);
U=get_umatrix_number_basis(H0,NT,F,Q,N);
[phin,En]=eig(U);
[phin,Es]=get_schur_ordered(N,En,phin,set_efn) ; %Reorder the schur vectors
mu=real(log(Es(1:num_efn)));
ds(itt)=abs(mu(1)-mu(2))
end
toc
figure(1)
hold on
plot(NTR,ds,'k.-','Markersize',5)
xlabel('NT')
ylabel('ds')

