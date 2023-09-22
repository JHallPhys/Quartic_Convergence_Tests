%==========================================================================
% Function to order the eigenvalues of the schur triangular matrix in theh
% Schur decomposition such that the imaginary parts are in either ascending
% or descending order. 
%   Inputs
%       * N (integer) size of the Schur matrices
%       * E_in(NxN array): upper triangular matrix with eigs on diagonal 
%       * Psi_in (NxN array): schur vectors corresoonding to E_in
%       * efn_set (string): String that indicates the two orderings
%   Output
%       * Psi_out (NxN array): Reordered Schur vectors spanning new
%       subspace
%       * E_out (1xN) array corresponding to the reordered spectrum
%==========================================================================


function [Psi_out,E_out] = get_schur_ordered(N,E_in,Psi_in,efn_set)


% Get order of the eigenvalues appearing on diagonal of E_in=exp(-1i*eps). 

lambda=ordeig(E_in); 

% Calculate the quasienergy spectrum E=1i*ln(lambda)

E=1i*log(lambda); 


%==========================================================================
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Note sort imag(E) needs to output the opposite order to the one you want!
% I can't remember where the documantation is for this. Cf with the order
% of the mu_n and the output E_out.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%==========================================================================
if isequal(efn_set,'G') % mu_1<mu_2<...
    
    [~,ind]=sort(imag(E)); 
    dummy(ind)=1:N;
    [Psi_out,E_out]=ordschur(Psi_in,E_in,dummy);
    E_out=diag(E_out);
    
end 
        

if isequal(efn_set,'L') %mu_1>mu_2>...
    
    [~,ind]=sort(imag(E),'descend'); 
    dummy(ind)=1:N;
    [Psi_out,E_out]=ordschur(Psi_in,E_in,dummy);
    E_out=diag(E_out);
    
end