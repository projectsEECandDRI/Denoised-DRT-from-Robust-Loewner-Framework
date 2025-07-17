function [Ek, Ak, Bk, Ck,s_LLs1] = state_space_mod(L,Ls,V,W,k)
% This function is part of the algorithm presented in the following
% publication:
% 
% Bansidhar Patel, Antonio Sorrentino, Ion Victor Gosea, 
% Athanasios C. Antoulas, and Tanja Vidakovic-Koch, "A data-driven, 
% noise-resilient algorithm for extraction of distribution of relaxation 
% times using the Loewner framework," Journal of Power Sources, vol. xxx, 
% Art. no. 237909, 2025. 
% Available: https://doi.org/10.1016/j.jpowsour.2025.237909

% % Singular value decomposition of the Loewner pencil (L,Ls)
[Y_LLs1,svd_LLs1, X_LLs1] = svd([L Ls],'econ');
s_LLs1 = diag(svd_LLs1);

[Y_LLs2,svd_LLs2, X_LLs2] = svd([L;Ls],'econ');
s_LLs2=diag(svd_LLs2);

Yk = Y_LLs1(:,1:k);
Xk = X_LLs2(:,1:k);

% % Singular value decomposition of the Loewner matrix
% [Y_L,svd_L,X_L]=svd(L);
% s_L=diag(svd_L);

% Yk = Y_L(:,1:k);
% Xk = X_L(:,1:k);

% % Reduced state space model interpolating the data

% % The reduced state-space model can be computed using either the singular
% % vectors (Xk and Yk) of the Loewner matrix L or the Loewner pencil (L,Ls).

Ek=-Yk.'*L*Xk;
Ak=-Yk.'*Ls*Xk;
Bk=Yk.'*V;
Ck=W.'*Xk;

end