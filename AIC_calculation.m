function AICc = AIC_calculation(R_i, tau_i, F, Z, k)
% This function is part of the algorithm presented in the following
% publication:
% 
% Bansidhar Patel, Antonio Sorrentino, Ion Victor Gosea, 
% Athanasios C. Antoulas, and Tanja Vidakovic-Koch, "A data-driven, 
% noise-resilient algorithm for extraction of distribution of relaxation 
% times using the Loewner framework," Journal of Power Sources, vol. xxx, 
% Art. no. 237909, 2025. 
% Available: https://doi.org/10.1016/j.jpowsour.2025.237909

% AIC_CALCULATION Computes the corrected Akaike Information Criterion (AICc)
%   for a series of models of order 1 to k using relaxation parameters.
%
% Inputs:
%   R_i   - [k x 1] Vector of relaxation resistances
%   tau_i - [k x 1] Vector of relaxation times
%   F     - [N x 1] Frequency vector
%   Z     - [N x 1] Measured complex impedance
%   k     - Scalar: maximum number of model terms
%
% Output:
%   AICc  - [k x 1] Corrected AIC values for model orders 1 to k

% Ensure column vectors
F = F(:); Z = Z(:); R_i = R_i(:); tau_i = tau_i(:);

% Number of frequency points
N = length(F);
omega = 2 * pi * F;  % Angular frequency

% Preallocate result
AICc = zeros(k, 1);

% Loop over model orders
for rr = 1:k
    % Select first rr parameters
    R_sub = R_i(1:rr);
    tau_sub = tau_i(1:rr);

    % Construct model response Hr using vectorized computation
    denom = 1 + 1i * omega * tau_sub.';      % [N x rr]
    Hr_matrix = R_sub.' ./ denom;           % [N x rr]
    Hr = sum(Hr_matrix, 2);                 % [N x 1] sum over rr terms

    % Compute residual sum of squares (RSS)
    residual = Z - Hr;
    RSS = sum(real(residual).^2 + imag(residual).^2);

    % Compute AIC and AICc
    AIC = 2 * rr + N * log(RSS);
    AICc(rr) = AIC + (2 * rr * (rr + 1)) / (N - rr - 1);
end

end
