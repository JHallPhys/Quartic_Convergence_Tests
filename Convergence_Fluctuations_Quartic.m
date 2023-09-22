clear all
close all

%==========================================================================
%   Initalise creation and anhllation Operators 
%==========================================================================

Nmin=100; % Max dimension for convergence 
Nmax=250; % Max dimension for convergence
dN=10;
Nrange=Nmin:dN:Nmax;
n_efn=20;
hbar_eff=0.2;
F=0.3
Qmax=1.5;
Pmax=1.5;
alpha=1;
gamma=0.001
beta=0;
set_efn='G'; % Invarient Subspace: Gain ('G') Loss ('L')
NT=100; % dt=2*pi/NT in integrator
% Rescaling
alpha=alpha*hbar_eff^2;
F=F/hbar_eff;
Qmax=Qmax/hbar_eff;
Pmax=Pmax/hbar_eff;
num_con=n_efn; % Number of eigenvalues for convergenvce
count=0;
for N=Nrange
   N
    count=count+1;
    [a,Q,P]=init_number_basis(N,1,1,1); % Get the operators in the number basis
    H0=0.5*(P^2)-0.5*beta*Q^2+0.25*alpha*(Q^4)-0.5*1i*gamma*P^2;
    U=get_umatrix_number_basis(H0,NT,F,Q,N);
    [phin,En]=eig(U);
    [phin,Es]=get_schur_ordered(N,En,phin,set_efn) ; %Reorder the schur vectors
    mu(count,:)=real(log(Es(1:num_con)));




end




mu_mean=NaN(length(Nrange),n_efn);
% Caluclate the N dependent mean and the varience
for j = 1:length(Nrange)
mu_mean(j,:)=sum(mu(1:j,:))./j;
mu_sigma(j,:)=abs(mu(j,:)-mu_mean(j,:));
end

mu_sigma=mu_sigma./gamma;
figure(1)
clf
hold on
plot(Nrange(2:end),mu_sigma(2:end,:),'k.-','Markersize',5)
for j =1:length(Nrange)-1
sigma_max(j)=max(mu_sigma(j+1,:));
end
plot(Nrange(2:end),sigma_max,'r.-','Markersize',10)
xlabel('N')
ylabel('\sigma \gamma^{-1}')