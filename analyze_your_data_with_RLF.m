% analyze_your_data_with_RLF.m
%
% This script provides an efficient and user-friendly pipeline to extract a 
% denoised Distribution of Relaxation Times (DRT) from your impedance data 
% using the Robust Loewner Framework (RLF) approach. 
% 
% This script implements the algorithm presented in the following 
% publication:
% 
% Bansidhar Patel, Antonio Sorrentino, Ion Victor Gosea, 
% Athanasios C. Antoulas, and Tanja Vidakovic-Koch, "A data-driven, 
% noise-resilient algorithm for extraction of distribution of relaxation 
% times using the Loewner framework," Journal of Power Sources, vol. xxx, 
% Art. no. 237909, 2025. 
% Available: https://doi.org/10.1016/j.jpowsour.2025.237909
% 
% In this implementation, the model order is selected based on the 
% ScreeNOT criterion, which enables optimal singular value thresholding 
% in the presence of correlated noise. 
% The MATLAB implementation of ScreeNOT used in this script is publicly 
% available at:
%
% Romanov, E. (2023). *ScreeNOT: MATLAB Implementation of the Scree Noise 
% Optimal Thresholding Method*. Available online: 
% https://github.com/eladromanov/ScreeNOT
% 
% --- Main Steps ---
% - Import Impedance Data
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
% - Reproduce the impedance data using RLF-derived Denoised DRT
% - Compute relative residuals for validation
% - Plot figures inclusing the impedance spectrum and DRT for comparison
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

%% Importing impedance data

% Load the frequency vector F [Hz] and complex impedance data Z

% Load your own EIS dataset or explore using one of the example datasets
% provided below

% load('sampleEISData_withNoise_from_R_2RC_ECM.mat'); 
load('sampleEISData_withNoise_from_R_2RC_RL_ECM.mat'); 
% load('sampleEISData_withNoise_from_R_RQ_RQ_ECM.mat'); 
% load('sampleEISData_withNoise_from_Randles_ECM.mat'); 

% Ensure frequency (F) and impedance (Z) data are column vectors
F = F(:); Z = Z(:);

%% Construction of the Loewner framework matrices and vectors 
% NOTE: Option `'real is set if one deals with data sets with imaginary
% entries
[L,Ls,V,W] = Loewner_Framework(1i*2*pi*F,Z,'real');

%% Determining the model order for SVD truncation 
% Estimate the optimal model order using ScreeNOT-based SV truncation.
% The imputation-based CDF strategy ('i') is used.
% The upper bound for the rank is set as: floor (min(size([L Ls])) / 2) - 1.
k_upper=floor(min(size([L Ls],1),size([L Ls],2))/2)-1; 
[Xest, Topt, k_LF] = ScreeNOT([L Ls], k_upper, 'i');

%% Construction of the reduced state-space model 
[Ek, Ak, Bk, Ck, s_LLs1] = state_space_mod(L,Ls,V,W,k_LF);

%% DRT calculation 
% Compute Distribution of Relaxation Times (DRT) from the model
[R_i_LF,tau_i_LF,pol_i_LF,res_i_LF] = DRT_values(Ek,Ak,Bk,Ck);

%% AIC-based filtering (Akaike information criterion) 

% Step 1: Sort DRT elements based on residue norms
% This helps identify the most influential poles (those with larger residues)
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

% Note: The final choice between r_min and r_drop may depend on the visual 
% inspection of AICc trends (e.g., see Figure 5).

% Final selection of optimal model order (number of relevant DRT elements)
% % Option 1: Based on the largest drop in successive AICc values
% r = r_drop;
% % Option 2: Based on the minimum AICc value
r = r_min;
% % Option 3: Based on visual inspection of AICc trends (e.g., see Figure 5)
% r = ; % Manually assigned from AICc trends

% Denoised DRT elements
R_i_RLF  = R_i_LF(1:r);
tau_i_RLF = tau_i_LF(1:r);

%% Reproduction of impedance data using RLF-derived DRT
% % Case 1: Reconstruct Impedance by Explicitly Separating R_ohm
% % (Assumes R_i_RLF(1) corresponds to R_ohm)

% Intrapolation
Hr_DRT=[];
    for i1=1:length(F)
        hr=0;
    for j1=2:r
    a= R_i_RLF(j1)./(1+1i*2*pi*F(i1)*tau_i_RLF(j1));
        hr=hr+a;
    end
        Hr_DRT(i1,1)=hr+R_i_RLF(1);
    end

%Extrapolation
F_BIG = logspace(-2,6,8001); 
Hr_DRT_F_BIG=[];
    for i1=1:length(F_BIG)
        hr=0;
    for j1=2:r
    a= R_i_RLF(j1)./(1+1i*2*pi*F_BIG(i1)*tau_i_RLF(j1));
        hr=hr+a;
    end
        Hr_DRT_F_BIG(i1,1)=hr+R_i_RLF(1);
    end

% % Case 2: Reconstruct Impedance Without Explicitly Separating R_ohm

% % Intrapolation
% Hr_DRT=[];
%     for i1=1:length(F)
%         hr=0;
%     for j1=1:r
%     a= R_i_RLF(j1)./(1+1i*2*pi*F(i1)*tau_i_RLF(j1));
%         hr=hr+a;
%     end
%         Hr_DRT(i1,1)=hr;
%     end

% %Extrapolation
% F_BIG = logspace(-2,6,8001); 
% Hr_DRT_F_BIG=[];
%     for i1=1:length(F_BIG)
%         hr=0;
%     for j1=1:r
%     a= R_i_RLF(j1)./(1+1i*2*pi*F_BIG(i1)*tau_i_RLF(j1));
%         hr=hr+a;
%     end
%         Hr_DRT_F_BIG(i1,1)=hr;
%     end    
    
%% Relative residuals

% Relative residuals between noisy data and reconstructed data 
ReZ_rel_DRT=((real(Z)-real(Hr_DRT))./(abs(Z)));
ImZ_rel_DRT=((imag(Z)-imag(Hr_DRT))./(abs(Z)));

% Relative difference between noise-free and noisy data
ReZ_rel_data=((real(Z)-real(Z))./(abs(Z)));
ImZ_rel_data=((imag(Z)-imag(Z))./(abs(Z)));

%% Plotting Figures

figure(1)
% Figure 1: Nyquist plot
plot(real(Z),-imag(Z),'o','color',[0 0 0],'LineWidth',3,...
    'markersize',7,'DisplayName','Data')
hold on
plot(real(Hr_DRT),-imag(Hr_DRT),'-',...
    'color',[0 1 1],'markersize',7,'linewidth',3,'DisplayName',...
    'RLF-Intrapolate')
plot(real(Hr_DRT_F_BIG),-imag(Hr_DRT_F_BIG),'-.',...
    'color',[0 0 1],'markersize',7,'linewidth',3,'DisplayName',...
    'RLF-Extrapolate')
grid on; axis equal;
xlabel('Re(Z) [\Omega]','FontSize',18,'FontWeight','normal')
ylabel('-Im(Z) [\Omega]','FontSize',18,'FontWeight','normal')
title('Nyquist plot','FontSize',18,'FontWeight','normal')
legend('Location','NorthWest','FontSize',14,'FontWeight','normal',...
    'NumColumns',1)

figure(2)
subplot(2,1,1)
% Figure 2A: DRT plot-loglog
stem(abs(tau_i_RLF),abs(R_i_RLF),...
    '-.p','filled','color',[0 0 1],'markersize',7,'LineWidth',2,...
    'DisplayName','RLF') 
set(gca,'xscal','log','yscal','log')
xlabel('\tau [s]','FontSize',18,'FontWeight','normal')
ylabel('|R_i| [\Omega]','FontSize',18,'FontWeight','normal')
title('DRT plot - loglog','FontSize',18,'FontWeight','normal')
grid on;
legend('Location','NorthWest','FontSize',14,'FontWeight','normal',...
    'NumColumns',1)

subplot(2,1,2)
% Figure 2B: DRT plot-semilogx
stem(abs(tau_i_RLF),real(R_i_RLF),...
    '-.p','filled','color',[0 0 1],'markersize',7,'LineWidth',2,...
    'DisplayName','RLF') 
set(gca,'xscal','log')
xlabel('\tau [s]','FontSize',18,'FontWeight','normal')
ylabel('R_i [\Omega]','FontSize',18,'FontWeight','normal')
title('DRT plot - semilogx','FontSize',18,'FontWeight','normal')
grid on;
legend('Location','NorthWest','FontSize',14,'FontWeight','normal',...
    'NumColumns',1)

figure(3)
% Figure 3: Relative residuals
plot(F,abs(ReZ_rel_DRT),':+','color',[0 0 1],...
    'markersize',7,'LineWidth',2,'DisplayName','Real-RLF')
hold on
plot(F,abs(ImZ_rel_DRT),'--o','color',[0 0 1],...
    'markersize',7,'LineWidth',2,'DisplayName','Imag-RLF')
set(gca,'xscal','log','yscal','log')
xlabel('Frequency [Hz]','FontSize',18,'FontWeight','normal')
ylabel('Residual [-]','FontSize',18,'FontWeight','normal')
title('Relative residuals','FontSize',18,'FontWeight','normal')
grid on;
ylim([1e-17 1e5]);
yticks([1e-15 1e-10 1e-5 1e0 1e5]);
yticklabels({'10^{-15}','10^{-10}','10^{-5}','10^{0}','10^{5}'});
legend('Location','NorthWest','FontSize',14,'FontWeight','normal',...
    'NumColumns',1)

figure(4)
% Figure 4: Singular value decay of the Loewner pencil ([L, Ls])
semilogy(s_LLs1./s_LLs1(1),':ok','MarkerSize',7,'linewidth',2,...
    'DisplayName','SVs of [L, Ls]')
hold on
semilogy(k_LF,s_LLs1(k_LF)/s_LLs1(1),'pr','MarkerSize',7,...
    'DisplayName','ScreeNOT (r_{opt})')
axis tight;grid on;
xlabel('Order of the model (r)')
ylabel('Singular value (SV)')
title('Singular value decay of the Loewner pencil','FontSize',18,...
    'FontWeight','normal')
legend('Location','NorthEast','FontSize',14,'FontWeight',...
    'normal','NumColumns',1)

figure(5)
% Figure 5: AIC-based filtering
plot(AICc,'^-r','MarkerSize',6,'linewidth',2,'DisplayName','AICc Values')
hold on
plot(r,AICc(r),'pb','MarkerSize',7,'DisplayName',...
    'Refined Optimum Order ({r_{opt}}*)')
xlabel('Order of the model (r)')
ylabel('AICc')
legend('Location','NorthEast','FontSize',14,'FontWeight',...
    'normal','NumColumns',1)

yyaxis right
plot(1.5:length(diff_AICc)+0.5,diff_AICc,'^-','MarkerSize',6,...
    'linewidth',2,'DisplayName','AICc Value Differences')
xlabel('Order of the model (r)')
ylabel('AICc differences')
legend('Location','NorthEast','FontSize',14,'FontWeight',...
    'normal','NumColumns',1)
title('AICc trends','FontSize',18,'FontWeight','normal')

