function [R_i,tau_i,pol,res] = DRT_values(E,A,B,C)
% This function is part of the algorithm presented in the following
% publication:
% 
% Bansidhar Patel, Antonio Sorrentino, Ion Victor Gosea, 
% Athanasios C. Antoulas, and Tanja Vidakovic-Koch, "A data-driven, 
% noise-resilient algorithm for extraction of distribution of relaxation 
% times using the Loewner framework," Journal of Power Sources, vol. xxx, 
% Art. no. 237909, 2025. 
% Available: https://doi.org/10.1016/j.jpowsour.2025.237909

% % Eigenvalue decomposition of 
% [U,T] = eig(A,E); % eigenvalue decomposition

% % Determination of the poles
% pol   = diag(T);     

% % Determination of the residues
% Bt = U\(E\B); 
% Ct = C*U;
% res = Bt.*Ct.';   

% % Determination resistances R_i and related time constants values tau_i 
% % constituting the DRT
% R_i   = (-res./pol);
% tau_i = abs(-1./pol);

% % One can also use the following eigenvalue decomposition for the direct  
% % determination of resistances R_i and related time constant values tau_i.  
% % This allows direct computation of resistances R_i and 
% % time constants tau_i, bypassing the need for poles and residues.  
% % It also circumvents the explicit computation of very large poles 
% % near zero time constants.
% % If no ohmic resistance (R0 = 0) is present in the impedance dataset,  
% % it is recommended to use following EVD to prevent numerical errors caused by  
% % computations at machine precision.  
% % Such numerical artifacts may arise due to approximations in the  
% % singular value decomposition (SVD) of the Loewner matrices, 
% % particularly when the impedance data approaches zero at high frequencies.  
% % This adjustment does not compromise the applicability of the method
% % or the accuracy of the analysis.


[UU,TT] = eig(E,-A);
Ctt = C*UU;    
Btt = UU\(-A\B); 
tau_i= diag(TT);
tau_i=tau_i(:);
R_i = Ctt.*Btt.';
R_i=R_i(:);

% % Determination of the poles and the residues
pol=-1./tau_i;
res=-R_i.*pol;

end
