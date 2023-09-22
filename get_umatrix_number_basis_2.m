%==========================================================================
% Function to calculate the Floquet matrix of a time-dependent,
% time-periodic Hamiltonian in the number basis.
% Input
%       * H0 (array NxN) : dimension of the husimi grid 
%       * NT (integer) : Number of time steps
%       * F (real) : Hamiltonian force parameter 
%       * Qop(NxN) : Position Operator
%   Output
%       * Uout (NxN array): Floquet matrix
%==========================================================================
function Uout=get_umatrix_number_basis(H0,NT,F,Q_op,N)

omega=1; % Reference Harmonic oscillator frequency
Omega=1; % Forcing period
hbar=1;

% reverseStr = ''; % String for counter
% display(' ') % Formating in terminal

T = 2*pi/omega; 
tstep = T/NT; 
Uout = eye(N); % U(0,0)=id
count=0; 
for  t=tstep:tstep:T
    count=count+1;
   
%     msg = sprintf('Constructing the Floquet matrix %d/%d', count, NT); % Tell you how many are done
%     fprintf([reverseStr, msg]);
%     reverseStr = repmat(sprintf('\b'), 1, length(msg));
    Uout = expm(-1i*(H0-Q_op*F*cos(Omega*t))*tstep/hbar) * Uout;
end 
 
end