% BoxPlot_100SynExps_R_RQ_RQ_ECM_RLF.m
%
% This script generates box plots showing the statistical distribution of
% impedance parameters extracted from the Robust Loewner Framework (RLF)-
% derived Distribution of Relaxation Times (DRT) for an R(RQ)(RQ) 
% Equivalent Circuit Model (ECM) containing two Zarc elements in series
% with an ohmic resistance.
%
% The analysis is performed across different noise levels, sampling
% densities, and frequency ranges. For each scenario, 100 synthetic 
% experiments are conducted.
% 
% This example script reproduces the results shown in Figures S7, S11, S12,
% and S16 of the paper:
%
% Bansidhar Patel, Antonio Sorrentino, Ion Victor Gosea, 
% Athanasios C. Antoulas, and Tanja Vidakovic-Koch, "A data-driven, 
% noise-resilient algorithm for extraction of distribution of relaxation 
% times using the Loewner framework," Journal of Power Sources, vol. xxx, 
% Art. no. 237909, 2025. 
% Available: https://doi.org/10.1016/j.jpowsour.2025.237909
%
% The steps implemented here correspond to the algorithm described in the 
% Supplementary Information of the above paper. In this implementation,
% the model order is selected based on the ScreeNOT criterion, which 
% enables optimal singular value thresholding in the presence of 
% correlated noise. 
% The MATLAB implementation of ScreeNOT used in this script is publicly 
% available at:
%
% Romanov, E. (2023). *ScreeNOT: MATLAB Implementation of the Scree Noise 
% Optimal Thresholding Method*. Available online: 
% https://github.com/eladromanov/ScreeNOT
%
% The ScreeNOT method is based on the following publication:
%
% Donoho, D., Gavish, M., and Romanov, E. (2023). 
% "ScreeNOT: Exact MSE-optimal singular value thresholding in correlated noise," 
% *The Annals of Statistics*, 51(1), 122–148.

% --- Main Steps ---
% - Generate synthetic (noise-free) impedance data from an R(RQ)(RQ) ECM
% - Add synthetic noise to the impedance data
% - Construct Loewner matrix (L), shifted Loewner matrix (Ls), and vectors 
%   (V and W) using 'Loewner_Framework.m'
% - Determine the model order using optimal SVD truncation via the ScreeNOT 
%   methodology, implemented through the 'ScreeNOT.m'.
% - Generate reduced-order state-space model (Ek, Ak, Bk, Ck) using
%   'state_space_mod.m'
% - Extract resistances R_i and time constants tau_i from state-space 
%   matrices using 'DRT_values.m'
% - AIC-based filtering (Akaike Information Criterion) using
%   'AIC_calculation.m' for refined optimal model order selection that 
%   identifies the most relevant DRT elements (denoised DRT)
% - Extract impedance parameters from RLF-derived Denoised DRT
% - Perform statistical analysis across synthetic experiments for 
%   parameter robustness
% - Visualize results using box plots for performance evaluation
% 
% Last revised: 16.07.2025
%-------------------------------------------------------------------------%
clear all; close all;

% Set default figure and plot properties for better visualization
set(0,'DefaultFigurePosition', [1,41,1280,607.33]);
set(0,'defaultlinelinewidth',3)
set(0,'defaultlinemarkersize',7)
set(0,'defaultaxesfontsize',18)
set(0,'DefaultAxesTitleFontWeight','normal');
set(0,'DefaultaxesFontWeight','normal') 
% set(0,'DefaultaxesFontWeight','bold') 
set(0,'DefaultTextFontName','Arial');
set(0,'defaultAxesXGrid','on') 
set(0,'defaultAxesYGrid','on')

%% Creating an impedance spectra data set of R(RQ)(RQ) ECM

% Define ECM parameters
R0=10;
R1=20;
tau1=1e-2;
phi1=0.85;
R2=150;
tau2=1;
phi2=0.85;

% Note: The sampling density (e.g., points per decade) and frequency range 
% can be modified to observe their effects on the performance of the RLF 
% approach in recovering ECM parameters.
% 
% Specifically, the following frequency vectors are used to reproduce 
% the corresponding figures: 
% Figure S7 with f = logspace(-2, 6, 81).', 
% Figure S11 with f = logspace(-2, 6, 41).', 
% Figure S12 with f = logspace(-2, 6, 25).',and 
% Figure S16 with f = logspace(-1,5,61).'.

f = logspace(-2,6,81).'; % frequency range of the impedance spectra [Hz]
w   = 2*pi*f; % Angular frequency

% R(RQ)(RQ) ECM impedance
Z_0 = R0+R1./(1+(1i*w*tau1).^phi1)+R2./(1+(1i*w*tau2).^phi2);

%% Initialize noise levels, experiment count, and preallocate result arrays

Noise_Level=[0 0.001, 0.005,0.01];  % Define noise levels (0 = noise-free)
n_exp=100;  % Number of synthetic experiments

R_i_all=nan(80,n_exp,length(Noise_Level));  
tau_i_all=nan(80,n_exp,length(Noise_Level));

R0_RLF_all = nan(n_exp,length(Noise_Level));
R1_RLF_all = nan(n_exp,length(Noise_Level));
tau1_RLF_all = nan(n_exp,length(Noise_Level));
R2_RLF_all = nan(n_exp,length(Noise_Level));
tau2_RLF_all = nan(n_exp,length(Noise_Level));

mean_values = nan(length(Noise_Level),5);
std_values = nan(length(Noise_Level),5);

%% Perform synthetic experiments across defined noise levels
for kk=1:length(Noise_Level)  
    rng('default'); % Set random seed for reproducibility
    for jj=1:n_exp

        %% Add synthetic noise to impedance data
        NL=Noise_Level(kk); 
        ReZ = real(Z_0) + (NL*abs(Z_0).*randn(length(Z_0),1));
        ImZ = imag(Z_0) + (NL*abs(Z_0).*randn(length(Z_0),1));
        Z_n=ReZ+1i.*ImZ; 
        
        % Store frequency and impedance data
        F = f(:); Z = Z_n(:);

%% Construction of the Loewner framework matrices and vectors 
        % NOTE: Option 'real' is set if one deals with data sets with 
        % imaginary entries
        [L,Ls,V,W] = Loewner_Framework(1i*2*pi*F,Z,'real');

        %% Determining the model order for SVD truncation 
        % Estimate the optimal model order using ScreeNOT-based 
        % SV truncation.
        % The imputation-based CDF strategy ('i') is used.
        % The upper bound for the rank is set as: 
        % floor (min(size([L Ls])) / 2 ) - 1.
        k_upper=floor(min(size([L Ls],1),size([L Ls],2))/2)-1; 
        [Xest, Topt, k_LF] = ScreeNOT([L Ls], k_upper, 'i');

        %% Construction of the reduced state-space model 
        [Ek, Ak, Bk, Ck, s_LLs1] = state_space_mod(L,Ls,V,W,k_LF);

        %% DRT calculation 
        % Compute Distribution of Relaxation Times (DRT) from the model
        [R_i_LF,tau_i_LF,pol_i_LF,res_i_LF] = DRT_values(Ek,Ak,Bk,Ck);
        
        %% AIC-based filtering (Akaike information criterion) 

        % Step 1: Sort DRT elements based on residue norms
        % This helps identify the most influential poles 
        % (those with larger residues)
        [Res_i, ind_i] = sort(abs(res_i_LF), 'descend');

        % Step 2: Reorder resistance and time constant values accordingly
        R_i_LF  = R_i_LF(ind_i);
        tau_i_LF = tau_i_LF(ind_i);

        % Step 3: Compute corrected Akaike Information Criterion (AICc)
        AICc = AIC_calculation(R_i_LF, tau_i_LF, F, Z, k_LF);

        % Step 4a: Select optimal model order based on minimum AICc value
        [~, r_min] = min(AICc);  % Order with minimum AICc

        % Step 4b: Alternative selection based on largest AICc drop
        diff_AICc = abs(AICc(1:end-1) - AICc(2:end));
        [~, idx_drop] = max(diff_AICc);
        r_drop = idx_drop + 1;   % Order with largest successive AICc drop

        % Note: The final choice between r_min and r_drop may depend on 
        % the visual inspection of AICc trends. 
        % In this example, r is selected based on the largest drop in 
        % successive AICc values across all synthetic experiments. 
        % However, for more accurate results, it is recommended to inspect 
        % the AICc trends individually for each simulation.

        % Final selection of optimal model order (number of relevant DRT
        % elements)
        r = r_min; % Selected based on minimum AICc values
                    
        % Denoised DRT elements
        R_i_RLF  = R_i_LF(1:r);
        tau_i_RLF = tau_i_LF(1:r);

        %% Impedance parameters calculation from computed discrete DRT
        % The process of computing impedance parameters from the 
        % obtained discrete DRT is discussed in Section S3.

        % Storing DRT parameters in descending order of time constants
        [~,ind_ii] = sort(abs(tau_i_RLF),'descend');
        R_i_RLF = R_i_RLF(ind_ii);
        tau_i_RLF = tau_i_RLF(ind_ii);

        % Remove duplicates from complex conjugate pairs (if any)
        rounded_column_R = round(abs(R_i_RLF), 4);
        rounded_column_tau = round(abs(tau_i_RLF), 4);
        [unique_values_R, ia_R] = unique(rounded_column_R, 'stable');    
        [unique_values_tau, ia_tau] = unique(rounded_column_tau, 'stable'); 
        R_i_RLF = R_i_RLF(ia_R);
        tau_i_RLF = tau_i_RLF(ia_tau);
        
        % Store cleaned DRT parameters across simulations
        R_i_all(1:length(unique_values_R), jj, kk) = R_i_RLF;
        tau_i_all(1:length(unique_values_tau), jj, kk) = tau_i_RLF;

        % Calculate impedance parameters from DRT
        % R0: Ohmic resistance (typically corresponds to the smallest tau)
        R0_RLF_all(jj,kk) = abs(R_i_RLF(end));
        % R1: High-frequency resistance (tau in range [1e-5, 5e-2] s)
        R1_RLF_all(jj,kk) = abs(sum(R_i_RLF(abs(tau_i_RLF)>=1e-5 &...
            abs(tau_i_RLF)<=5e-2)));
        % tau1: Weighted average relaxation time for R1
        weights1=abs(R_i_RLF(abs(tau_i_RLF)>=1e-5 & abs(tau_i_RLF)...
            <=5e-2))./R1_RLF_all(jj,kk);
        tau1_RLF_all(jj,kk) = sum(abs(tau_i_RLF(abs(tau_i_RLF)>=1e-5...
            & abs(tau_i_RLF)<=5e-2)).*weights1)./ sum(weights1);
        % R2: Low-frequency relaxation resistance (tau > 5e-2 s)
        R2_RLF_all(jj,kk) = abs(sum(R_i_RLF(abs(tau_i_RLF)>5e-2)));
        % tau2: Weighted average relaxation time for R2
        weights2=abs(R_i_RLF(abs(tau_i_RLF)>5e-2))./R2_RLF_all(jj,kk);                  
        tau2_RLF_all(jj,kk) = sum(abs(tau_i_RLF(abs(tau_i_RLF)...
            >5e-2)).*weights2)./ sum(weights2);
                                                  
    end
        % Store impedance parameters for all noise levels
        impedance_parameters=[R1_RLF_all(:,kk),tau1_RLF_all(:,kk),...
            R2_RLF_all(:,kk),tau2_RLF_all(:,kk),R0_RLF_all(:,kk)];
   
        %% Statistical Parameter Calculation
        % Compute 25th (Q1) and 75th (Q3) percentiles for each parameter
        q1_values = prctile(impedance_parameters, 25, 1);
        q3_values = prctile(impedance_parameters, 75, 1);
        iqr_values = q3_values - q1_values;  % Interquartile range (IQR)
        
        % Define bounds for outlier detection
        lower_bound = q1_values - 1.5 * iqr_values;  
        upper_bound = q3_values + 1.5 * iqr_values;  
        
        % Logical matrix indicating outlier values
        outliers = (impedance_parameters < lower_bound) |...
            (impedance_parameters > upper_bound);
        
        % Remove outliers by assigning NaN
        no_outliers = impedance_parameters;
        no_outliers(outliers) = NaN;
        
        % Compute standard deviation and mean (ignoring NaNs)
        std_values(kk,:)= std(no_outliers, 0, 1, 'omitnan');
        mean_values(kk,:) = mean(no_outliers, 1, 'omitnan');
end

% True impedance/ECM parameter values for comparison
true_values = [R1,tau1,R2,tau2,R0];

%% Box Plots for All Noise Levels
figure(1);
% Figure 1A: Box plot of R1
subplot('Position',[0.079,0.53,0.4,0.4]);
boxplot(R1_RLF_all(:,2:end), ...
    'Labels', arrayfun(@num2str, Noise_Level(2:end), 'UniformOutput',...
    false),'Notch', 'on', 'Symbol', 'r+');
hold on
plot(1:length(Noise_Level(2:end)), mean_values(2:end,1), 'co', ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Mean');
errorbar(1:length(Noise_Level(2:end)), mean_values(2:end,1),...
    std_values(2:end,1),'c', 'LineStyle', 'none', 'CapSize', 10,...
    'LineWidth', 2, 'DisplayName', 'Std. Dev.');
yline(true_values(1), 'k--', 'LineWidth', 3, 'DisplayName', 'True Value');
ylabel('R_1 [\Omega]','FontSize',18,'FontWeight','normal');
xlabel('Noise Level (NL)','FontSize',18,'FontWeight','normal');
ylim([10,30])
grid on;
legend('Location','NorthWest','FontSize',14,'FontWeight','normal',...
    'NumColumns',3)
annotation('textbox','String',{'A'},'FontSize',18,'FontWeight','bold',...
    'LineStyle','none','Position',[0.455,0.87,0.025,0.06]);

% Figure 1B: Box plot of tau1
subplot('Position',[0.573,0.53,0.4,0.4]);
boxplot(tau1_RLF_all(:,2:end), ...
    'Labels', arrayfun(@num2str, Noise_Level(2:end), 'UniformOutput',...
    false),'Notch', 'on', 'Symbol', 'r+');
hold on
plot(1:length(Noise_Level(2:end)), mean_values(2:end,2), 'co', ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Mean');
errorbar(1:length(Noise_Level(2:end)), mean_values(2:end,2),...
    std_values(2:end,2),'c', 'LineStyle', 'none', 'CapSize', 10,...
    'LineWidth', 2, 'DisplayName', 'Std. Dev.');
yline(true_values(2), 'k--', 'LineWidth', 3, 'DisplayName', 'True Value');
ylabel('\tau_1 [s]','FontSize',18,'FontWeight','normal');
xlabel('Noise Level (NL)','FontSize',18,'FontWeight','normal');
ylim([0.001,0.03])
grid on;
legend('Location','NorthWest','FontSize',14,'FontWeight','normal',...
    'NumColumns',3)
annotation('textbox','String',{'B'},'FontSize',18,'FontWeight','bold',...
    'LineStyle','none','Position',[0.95,0.87,0.025,0.06]);

figure(2);
% Figure 2A: Box plot of R2 
subplot('Position',[0.079,0.53,0.4,0.4]);
boxplot(R2_RLF_all(:,2:end), ...
    'Labels', arrayfun(@num2str, Noise_Level(2:end), 'UniformOutput',...
    false),'Notch', 'on', 'Symbol', 'r+');
hold on
plot(1:length(Noise_Level(2:end)), mean_values(2:end,3), 'co', ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Mean');
errorbar(1:length(Noise_Level(2:end)), mean_values(2:end,3),...
    std_values(2:end,3), 'c', 'LineStyle', 'none', 'CapSize', 10,...
    'LineWidth', 2, 'DisplayName', 'Std. Dev.');
yline(true_values(3), 'k--', 'LineWidth', 3, 'DisplayName', 'True Value');
ylabel('R_2 [\Omega]','FontSize',18,'FontWeight','normal');
xlabel('Noise Level (NL)','FontSize',18,'FontWeight','normal');
ylim([120,170])
grid on;
legend('Location','NorthWest','FontSize',14,'FontWeight','normal',...
    'NumColumns',3)
annotation('textbox','String',{'C'},'FontSize',18,'FontWeight','bold',...
    'LineStyle','none','Position',[0.455,0.87,0.025,0.06]);

% Figure 2B: Box plot of tau2 
subplot('Position',[0.573,0.53,0.4,0.4]);
boxplot(tau2_RLF_all(:,2:end), ...
    'Labels', arrayfun(@num2str, Noise_Level(2:end), 'UniformOutput',...
    false),'Notch', 'on', 'Symbol', 'r+');
hold on
plot(1:length(Noise_Level(2:end)), mean_values(2:end,4), 'co', ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Mean');
errorbar(1:length(Noise_Level(2:end)), mean_values(2:end,4),...
    std_values(2:end,4),'c', 'LineStyle', 'none', 'CapSize', 10,...
    'LineWidth', 2, 'DisplayName', 'Std. Dev.');
yline(true_values(4), 'k--', 'LineWidth', 3, 'DisplayName', 'True Value');
ylabel('\tau_2 [s]','FontSize',18,'FontWeight','normal');
xlabel('Noise Level (NL)','FontSize',18,'FontWeight','normal');
ylim([0.1,3])
grid on;
legend('Location','NorthWest','FontSize',14,'FontWeight','normal',...
    'NumColumns',3)
annotation('textbox','String',{'D'},'FontSize',18,'FontWeight','bold',...
    'LineStyle','none','Position',[0.95,0.87,0.025,0.06]);

figure(3);
% Figure 3A: Box plot of R0 
subplot('Position',[0.079,0.53,0.4,0.4]);
boxplot(R0_RLF_all(:,2:end), ...
    'Labels', arrayfun(@num2str, Noise_Level(2:end), 'UniformOutput',...
    false), 'Notch', 'on', 'Symbol', 'r+');
hold on
plot(1:length(Noise_Level(2:end)), mean_values(2:end,5), 'co', ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Mean');
errorbar(1:length(Noise_Level(2:end)), mean_values(2:end,5),...
    std_values(2:end,5),'c', 'LineStyle', 'none', 'CapSize', 10,...
    'LineWidth', 2, 'DisplayName', 'Std. Dev.');
yline(true_values(5), 'k--', 'LineWidth', 3, 'DisplayName', 'True Value');
ylabel('R_0 [\Omega]','FontSize',18,'FontWeight','normal');
xlabel('Noise Level (NL)','FontSize',18,'FontWeight','normal');
ylim([9,11])
grid on;
legend('Location','NorthWest','FontSize',14,'FontWeight','normal',...
    'NumColumns',3)
annotation('textbox','String',{'E'},'FontSize',18,'FontWeight','bold',...
    'LineStyle','none','Position',[0.455,0.87,0.025,0.06]);

