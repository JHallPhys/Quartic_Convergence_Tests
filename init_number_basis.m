%==========================================================================
%  Function to write the Harmonic oscillator ladder operators and the
%  position and momentum operators in the number basis. Assuimg that
%  hbar=omega=1;
%  
%   Input 
%       * N (integer): Hilbert space dimension
%       * hbar,m,omega (reals): Natural units of the system
%   Output
%       * a_out (NxN array): Anhillation operator
%       * q_out (NxN array): Position operator
%       * p_out (NxN array): Momentum operator
%==========================================================================

function [a_out,q_out,p_out]=init_number_basis(N,hbar,m,omega)

    % Define Anhillation operator
    a_out = diag(sqrt(1:N-1),1); 
    % Define Position operator
    q_out=sqrt(hbar/(2*m*omega))*(a_out+a_out'); 
    % Define Momentum operator
    p_out=-1i*sqrt(0.5*hbar*m*omega)*(a_out-a_out');
    

end