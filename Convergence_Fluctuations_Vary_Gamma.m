clear all
close all

%==========================================================================
%   Initalise creation and anhllation Operators 
%==========================================================================
gamma_range=[0.001,0.01,0.05,0.1]
 tic
for itt_gamma=1:length(gamma_range)

Nmin=100; % Max dimension for convergence 
Nmax=180; % Max dimension for convergence
dN=5;
Nrange=Nmin:dN:Nmax;
n_efn=20;
hbar_eff=0.14;
F=0.3;
Qmax=1.5;
Pmax=1.5;
alpha=1;
beta=0;
gamma=gamma_range(itt_gamma);
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
   
    count=count+1;
    [a,Q,P]=init_number_basis(N,1,1,1); % Get the operators in the number basis
    H0=0.5*(P^2)-0.5*beta*Q^2+0.25*alpha*(Q^4)-0.5*1i*gamma*P^2;
    U=get_umatrix_number_basis(H0,NT,F,Q,N);
    [phin,En]=eig(U);
    [phin,Es]=get_schur_ordered(N,En,phin,set_efn) ; %Reorder the schur vectors
    mu(count,:)=real(log(Es(1:num_con)));




end



% figure(1)
% hold on
% for j = 1:num_con
% plot(Nrange,mu(:,j),'k.-','Markersize',5)
% end
% xlabel('N')
% ylabel('\mu')


mu_mean=NaN(length(Nrange),n_efn);
% Caluclate the N dependent mean and the varience
for j = 1:length(Nrange)
mu_mean(j,:)=sum(mu(1:j,:))./j;
mu_sigma(j,:)=abs(mu(j,:)-mu_mean(j,:));
end


% figure(1)
% clf
% hold on
% plot(Nrange(2:end),mu_sigma(2:end,:),'k.-','Markersize',5)
% for j =1:length(Nrange)-1
% sigma_max(j)=max(mu_sigma(j+1,:));
% end
% plot(Nrange(2:end),sigma_max,'r.-','Markersize',10)
% xlabel('N')
% ylabel('\sigma')


mu_sigma=mu_sigma./gamma;
for j =1:length(Nrange)-1
sigma_max(j,itt_gamma)=max(mu_sigma(j+1,:));
sigma_av(j,itt_gamma)=mean(mu_sigma(j+1,:));
end
toc
end



figure(1)
clf
hold on
plot(Nrange(2:end),sigma_max(:,1),'r.-','Markersize',10)
plot(Nrange(2:end),sigma_max(:,2),'b.-','Markersize',10)
plot(Nrange(2:end),sigma_max(:,3),'g.-','Markersize',10)
plot(Nrange(2:end),sigma_max(:,4),'m.-','Markersize',10)
xlabel('N')
ylabel('\sigma \gamma^{-1}')
legend('0.001','0.01','0.05','0.1')


figure(2)
clf
hold on
plot(Nrange(2:end),sigma_av(:,1),'r.-','Markersize',10)
plot(Nrange(2:end),sigma_av(:,2),'b.-','Markersize',10)
plot(Nrange(2:end),sigma_av(:,3),'g.-','Markersize',10)
plot(Nrange(2:end),sigma_av(:,4),'m.-','Markersize',10)
xlabel('N')
ylabel('<\sigma \gamma^{-1}>')
legend('0.001','0.01','0.05','0.1')
