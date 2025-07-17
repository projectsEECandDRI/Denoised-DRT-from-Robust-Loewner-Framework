function [L,Ls,V,W] = Loewner_Framework(w,H,varargin)
% This function is part of the algorithm presented in the following
% publication:
% 
% Bansidhar Patel, Antonio Sorrentino, Ion Victor Gosea, 
% Athanasios C. Antoulas, and Tanja Vidakovic-Koch, "A data-driven, 
% noise-resilient algorithm for extraction of distribution of relaxation 
% times using the Loewner framework," Journal of Power Sources, vol. xxx, 
% Art. no. 237909, 2025. 
% Available: https://doi.org/10.1016/j.jpowsour.2025.237909

% making it sure that the inputs are in a column vector
w = w(:); H = H(:);

N = length(w);

% check for the 'real' option
% NOTE: Option `'real is set if one deals with data sets with imaginary
% numbers.

option = 0;
if nargin == 3
   if strcmp(varargin{1},'real')
      option = 1;
   end
end

% ensuring that the inputs have an even number of elements for 
% constructing the model having real entries 
if option == 1
    if rem(N,2)~=0
        N=N-1;
    end
end

% alternate partition of the given measurements
la = w(1:2:N); W = H(1:2:N);
mu = w(2:2:N); V = H(2:2:N);

% adding the complex conjugate values to the left and right 
% subsets for the real model
if option == 1
    la=[la.';la'];
    W=[W.';W'];       
    mu=[mu.';mu'];    
    V=[V.';V'];       
    
    la=la(:); W=W(:); mu=mu(:); V=V(:);
end

% constructing the Loewner and shifted-Loewner matrices
for i=1:length(mu)
    for j=1:length(la)
        L(i,j)=(V(i)-W(j))/(mu(i)-la(j));
        Ls(i,j)=(mu(i)*V(i)-la(j)*W(j))/(mu(i)-la(j));
    end
end

% transforming the complex Loewner and shifted-Loewner matrices 
% to obtain matrices with real entries (See step 4 in supplementary
% information). 
if option == 1
    P=(1/sqrt(2))*[1 1i;1 -1i];
    PP=P;
    for i=1:N/2-1
        PP = blkdiag(PP,P);
    end  
    P=PP;

    L=P'*L*P;
    Ls=P'*Ls*P;
    W=(W.'*P).';
    V=P'*V;
    
% discarding possibly very small imaginary parts (due to round off)
  L=real(L); Ls=real(Ls); W=real(W); V=real(V);
  
end

end
