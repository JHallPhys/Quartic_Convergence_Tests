clear all
close all

%==========================================================================
%   Initalise creation and anhllation Operators 
%==========================================================================

Nmin=126; % Max dimension for convergence 
Nmax=130; % Max dimension for convergence
n_efn=3;
hbar_eff=0.2;
dN=1
F=0.3
Qmax=1.5;
Pmax=1.5;
alpha=1;
gamma=0.001
beta=0;
set_efn='G'; % Invarient Subspace: Gain ('G') Loss ('L')
NT=1000; % dt=2*pi/NT in integrator
% Rescaling
alpha=alpha*hbar_eff^2;
F=F/hbar_eff;
Qmax=Qmax/hbar_eff;
Pmax=Pmax/hbar_eff;


num_con=n_efn; % Number of eigenvalues for convergenvce
count=0;
for N=Nmin:dN:Nmax
   N
    count=count+1;
    [a,Q,P]=init_number_basis(N,1,1,1); % Get the operators in the number basis
    H0=0.5*(P^2)-0.5*beta*Q^2+0.25*alpha*(Q^4)-0.5*1i*gamma*P^2;
    U=get_umatrix_number_basis(H0,NT,F,Q,N);
    [phin,En]=eig(U);
    [phin,Es]=get_schur_ordered(N,En,phin,set_efn) ; %Reorder the schur vectors
    Cn(count,:)=real(log(Es(1:num_con)));

end


figure(1)
hold on
for j = 1:num_con
plot(Nmin:dN:Nmax,Cn(:,j),'k.-','Markersize',5)
end
xlabel('N')
ylabel('\mu')

