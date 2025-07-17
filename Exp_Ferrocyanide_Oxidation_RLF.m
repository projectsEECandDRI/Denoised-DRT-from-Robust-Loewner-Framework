% Exp_Ferrocyanide_Oxidation_RLF.m
%
% Extraction of the tempered Distribution of Relaxation Times (DRT) from 
% experimental impedance data of ferrocyanide oxidation at different 
% electrode rotation rates.
% 
% This example script reproduces the results shown in Figures 8C and 8D,  
% as well as Table S4 of the paper:
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
% 
% The imepdance data used in this script is publicly available at:
% 
% Patel, Bansidhar; Sorrentino, Antonio; Vidakovic-Koch, Tanja, 2025,
% "Data-Driven Analysis of Electrochemical Impedance Spectroscopy Using
% the Loewner Framework", https://doi.org/10.17617/3.MMLZAK, Edmond, V2 
%
% This data is part of the following publication:
% 
% B. Patel, A. Sorrentino, and T. Vidakovic-Koch, "Data-driven analysis of
% electrochemicalimpedance spectroscopy using the Loewner framework,"
% iScience, vol. 28, no. 3, p. 111987,2025.
% 
% --- Main Steps ---
% - Importing experimental impedance spectra data sets
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
% - Extract impedance parameters from RLF-derived Denoised DRT
% - Plot the impedance and DRT figures for comparison
% - Display extracted impedance parameters versus rotation rate
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
newcolors = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];...
    [0.9290 0.6940 0.1250]];

%% Importing experimental impedance spectra data sets 

EIS1=dlmread('experimentalEISData_Ferrocyanide_Oxidation_600rpm.dat', '\t');
EIS2=dlmread('experimentalEISData_Ferrocyanide_Oxidation_900rpm.dat', '\t');
EIS3=dlmread('experimentalEISData_Ferrocyanide_Oxidation_3200rpm.dat', '\t');

F1 = EIS1(:,1)  ; % frequency [Hz]
F2 = EIS2(:,1);
F3 = EIS3(:,1);

Z1=EIS1(:,2)+1j*-EIS1(:,3); % impedance [Ohm]
Z2=EIS2(:,2)+1j*-EIS2(:,3);
Z3=EIS3(:,2)+1j*-EIS3(:,3);

F_all={F1(:),F2(:),F3(:)};
Z_all={Z1(:),Z2(:),Z3(:)};

S=1.96e-5;       % surface are of the electrode [m2]
omega_R=[600;900;3200]; % RPM

%% Extraction of impedance parameters and DRTs across operating conditions
for i=1:3
 
F=eval(['F',num2str(i)]); % Frequency in Hz
Z=eval(['Z',num2str(i)]); % Impedance in Ohm
% Z=Z*S; % Impedance in Ohm m2

F=F(2:end);Z=Z(2:end);
%% Construction of the Loewner framework matrices and vectors 
% NOTE: Option `'real is set if one deals with data sets with imaginary
% entries
[L,Ls,V,W] = Loewner_Framework(1i*2*pi*F,Z,'real');

%% Determining the model order for SVD truncation 
% Estimate the optimal model order using ScreeNOT-based SV truncation.
% The imputation-based CDF strategy ('i') is used.
% The upper bound for the rank is set as: floor(min(size([L Ls]))/2) - 1.
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
% inspection of AICc trends (e.g., see Figure 4).

% Final selection of optimal model order (number of relevant DRT elements)
r = r_drop; % Selected based on the largest drop in successive AICc values

% Denoised DRT elements
R_i_RLF  = R_i_LF(1:r);
tau_i_RLF = tau_i_LF(1:r);

%% Reproduction of impedance data using RLF-derived DRT

Hr_DRT=[];
    for i1=1:length(F)
        hr=0;
    for j1=2:r
    a= R_i_RLF(j1)./(1+1i*2*pi*F(i1)*tau_i_RLF(j1));
        hr=hr+a;
    end
        Hr_DRT(i1,1)=hr+R_i_RLF(1);
    end

%% Relative residuals

% Relative residuals between noisy data and reconstructed data 
ReZ_rel_DRT=((real(Z)-real(Hr_DRT))./(abs(Z)));
ImZ_rel_DRT=((imag(Z)-imag(Hr_DRT))./(abs(Z)));

%% Impedance parameters calculation from computed discrete DRT
% The process of computing impedance parameters from the 
% obtained discrete DRT is discussed in Section S3.

% Storing DRT parameters in descending order of time constants
[~,ind_ii] = sort(abs(tau_i_RLF),'descend');
R_i_RLF = R_i_RLF(ind_ii);
tau_i_RLF = tau_i_RLF(ind_ii);

% Remove duplicates from complex conjugate pairs (if any)
% rounded_column_R = round(abs(R_i_RLF), 4);
% rounded_column_tau = round(abs(tau_i_RLF), 4);
% [unique_values_R, ia_R] = unique(rounded_column_R, 'stable');    
% [unique_values_tau, ia_tau] = unique(rounded_column_tau, 'stable'); 
% R_i_RLF = R_i_RLF(ia_R);
% tau_i_RLF = tau_i_RLF(ia_tau);

% Calculate impedance parameters from DRT
% R0: Ohmic resistance (typically represents high-frequency intercept)
R0_RLF(i,1)=abs(R_i_RLF(end));
% Z0: Resistance associated with finite-length Warburg (FLW) element 
% – sum of R_i at >= 5e-5 s (mass transport region)
Z0_RLF(i,1)=abs(sum(R_i_RLF(abs(tau_i_RLF)>=5e-5)));
% tau0: Characteristic time constant of FLW – approximated from maximum tau 
% in the FLW region, adjusted by the peak position factor (0.4053)
tau0_RLF(i,1)=max(abs(tau_i_RLF(abs(tau_i_RLF)>=5e-5)))/0.4053;
% R_ct: Charge transfer resistance – sum of R_i in the kinetic region 
% (1e-6 s >= & < 5e-5 s)
R_RLF(i,1)=abs(sum(R_i_RLF(abs(tau_i_RLF)>=1e-6 & abs(tau_i_RLF)<5e-5)));
% C_dl: Double-layer capacitance estimate using tau and R in the kinetic region
C_RLF(i,1)=max(abs(tau_i_RLF(abs(tau_i_RLF)>=1e-6 & abs(tau_i_RLF)<...
    5e-5)))./abs(sum(R_i_RLF(abs(tau_i_RLF)>=1e-6 & abs(tau_i_RLF)<5e-5)));     
  
%% Plotting Figures

figure(1)
% Figure 1A: Nyquist plot
subplot('Position',[0.065 0.55 0.40 0.40])
plot(real(Z),-imag(Z),'o','markersize',7,'DisplayName',...
    [num2str(omega_R(i)),' rpm'],'color',newcolors(i,:))
hold on;
plot(real(Hr_DRT),-imag(Hr_DRT),'-.b','linewidth',3,...
    'DisplayName','RLF')
xlabel('Re(Z) [\Omega]','FontSize',18,'FontWeight','normal')
ylabel('-Im(Z) [\Omega]','FontSize',18,'FontWeight','normal')
axis equal; grid on;
xlim([0 48]);
ylim([0 20])
xticks([0 10 20 30 40]);xticklabels({'0','10','20','30','40'});
yticks([0 10 20]);yticklabels({'0','10','20'});
legend('Location','NorthWest','FontSize',14,'FontWeight',...
    'normal','NumColumns',3)
annotation('textbox','String',{'C'},'FontSize',18,'FontWeight','bold',...
    'LineStyle','none','Position',[0.44,0.87,0.025,0.06]);

% Figure 1B: Distribution of relaxation times
subplot('Position',[0.57 0.575 0.4 0.35])
stem(abs(tau_i_LF(1:r)),abs(R_i_LF(1:r)),'-.p','filled',...
    'markersize',7,'LineWidth',2,'DisplayName',...
    [num2str(omega_R(i)),' rpm'],'color',newcolors(i,:)) 
hold on
set(gca,'xscal','log','yscal','log')
xlabel('\tau [s]','FontSize',18,'FontWeight','normal')
ylabel('\midR_{i}\mid [\Omega]','FontSize',18,'FontWeight','normal')
grid on;
xlim([1e-8 1e0]); ylim([1e-1 1e2]);
xticks([1e-8 1e-6 1e-4 1e-2 1e0]);
xticklabels({'10^{-8}','10^{-6}','10^{-4}','10^{-2}','10^{0}'});
yticks([1e-1 1e0 1e1 1e2]);
yticklabels({'10^{-1}','10^{0}','10^{1}','10^{2}'});
legend('Location','NorthWest','FontSize',14,'FontWeight','normal','NumColumns',3)
annotation('textbox','String',{'D'},'FontSize',18,'FontWeight','bold',...
    'LineStyle','none','Position',[0.945,0.87,0.025,0.06]);

figure(2)
% Figure 2: Relative residuals
plot(F,abs(ReZ_rel_DRT),':+','color',[0 0 1],...
    'markersize',7,'LineWidth',2,'DisplayName',['Real-RLF for ',...
    num2str(omega_R(i)),'rpm'],'color',newcolors(i,:))
hold on
plot(F,abs(ImZ_rel_DRT),'--o','color',[0 0 1],...
    'markersize',7,'LineWidth',2,'DisplayName',['Imag-RLF for ',...
    num2str(omega_R(i)),'rpm'],'color',newcolors(i,:))
set(gca,'xscal','log','yscal','log')
xlabel('Frequency [Hz]','FontSize',18,'FontWeight','normal')
ylabel('Residual [-]','FontSize',18,'FontWeight','normal')
grid on;
xlim([1e-2 1e5]); ylim([1e-8 1e2]);
xticks([1e-2 1e-1 1e0 1e1 1e2 1e3 1e4 1e5]);
xticklabels({'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}',...
    '10^{3}','10^{4}','10^{5}'});
yticks([1e-8 1e-6 1e-4 1e-2 1e0 1e2]);
yticklabels({'10^{-8}','10^{-6}','10^{-4}','10^{-2}','10^{0}','10^{2}'});
legend('Location','NorthWest','FontSize',14,'FontWeight','normal',...
    'NumColumns',3)

figure(3)
% Figure 3: Singular value decay of the Loewner pencil ([L, Ls])
semilogy(s_LLs1./s_LLs1(1),':o','MarkerSize',7,'linewidth',1,...
    'DisplayName',['SVs of [L, Ls] for ', num2str(omega_R(i)),'rpm'],...
    'color',newcolors(i,:))
hold on
semilogy(k_LF,s_LLs1(k_LF)/s_LLs1(1),'p','MarkerSize',8,...
    'DisplayName','ScreeNOT (r_{opt})','color',newcolors(i,:))
axis tight;grid on;
xlabel('Order of the model (r)')
ylabel('Singular value (SV)')
legend('Location','NorthEast','FontSize',14,'FontWeight',...
    'normal','NumColumns',1)

figure(4)
% Figure 4: AIC-based filtering
subplot(2,1,1)
plot(AICc,'^-','MarkerSize',6,'linewidth',1,'DisplayName',...
    ['AICc Values for ', num2str(omega_R(i)),'rpm'],...
    'color',newcolors(i,:))
hold on
plot(r,AICc(r),'p','MarkerSize',8,'DisplayName',...
    'Refined Optimum Order ({r_{opt}}*)','color',newcolors(i,:))
xlabel('Order of the model (r)')
ylabel('AICc')
legend('Location','NorthEast','FontSize',14,'FontWeight',...
    'normal','NumColumns',2,'Orientation','horizontal')

subplot(2,1,2)
plot(1.5:length(diff_AICc)+0.5,diff_AICc,'^-','MarkerSize',6,...
    'linewidth',1,'DisplayName',['AICc Value Differences for ',...
    num2str(omega_R(i)),'rpm'],'color',newcolors(i,:))
hold on
xlabel('Order of the model (r)')
ylabel('AICc differences')
legend('Location','NorthEast','FontSize',14,'FontWeight',...
    'normal','NumColumns',1)

end

%% Display table of extracted impedance parameters versus rotation rate
% Create a table with impedance parameters extracted at each rotation rate
table = table(omega_R, R0_RLF, R_RLF, C_RLF, Z0_RLF, tau0_RLF,...
    'VariableNames', {'Rotation speed [rpm]', 'R0', 'R_ct',...
    'C','R_d', 'tau_d'});
                    
% Display the table in the Command Window
disp(table);

