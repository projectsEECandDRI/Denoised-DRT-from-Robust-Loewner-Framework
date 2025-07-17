% Example_R_RQ_RQ_ECM_LF.m
%
% Extraction of the tempered Distribution of Relaxation Times (DRT) from 
% noisy impedance spectra using the Loewner Framework (LF), applied to
% an R(RQ)(RQ) Equivalent Circuit Model (ECM) containing two Zarc elements  
% in series with an ohmic resistance.
%
% This example script reproduces the results shown in Figure 4 and 
% Figure S3 of the paper:
%
% Bansidhar Patel, Antonio Sorrentino, Ion Victor Gosea, 
% Athanasios C. Antoulas, and Tanja Vidakovic-Koch, "A data-driven, 
% noise-resilient algorithm for extraction of distribution of relaxation 
% times using the Loewner framework," Journal of Power Sources, vol. xxx, 
% Art. no. 237909, 2025. 
% Available: https://doi.org/10.1016/j.jpowsour.2025.237909
%
% The steps implemented here correspond to the algorithm described in the 
% Supplementary Information of the above paper. Since the model order is 
% chosen based on the rank of the Loewner pencil, this implementation is 
% also consistent with the approach described in the earlier publication:
%
% Sorrentino, A., Patel, B., Gosea, I.V., Antoulas, A.C., and 
% Vidakovic-Koch, T. (2023). "Determination of the distribution of 
% relaxation times through Loewner framework: A direct and versatile 
% approach." *Journal of Power Sources*, 585, 233575.  
% https://doi.org/10.1016/j.jpowsour.2023.233575

% --- Main Steps ---
% - Generate synthetic (noise-free) impedance data from an R(RQ)(RQ) ECM
% - Add synthetic noise to the impedance data
% - Construct Loewner matrix (L), shifted Loewner matrix (Ls), and vectors 
%   (V and W) using 'Loewner_Framework.m'
% - Determine the model order via SVD truncation based on the 
%   rank of Loewner pencil ([L Ls])
% - Generate reduced-order state-space model (Ek, Ak, Bk, Ck) using
%   'state_space_mod.m'
% - Extract resistances R_i and time constants tau_i from state-space 
%   matrices using 'DRT_values.m'
% - Reproduce the impedance data using LF-derived DRT
% - Compute relative residuals for validation
% - Plot the impedance and DRT figures for comparison
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

f = logspace(-2,6,81).'; % frequency range of the impedance spectra [Hz]
w   = 2*pi*f; % Angular frequency

% R(RQ)(RQ) ECM impedance
Z_0 = R0+R1./(1+(1i*w*tau1).^phi1)+R2./(1+(1i*w*tau2).^phi2); 

%% Add synthetic noise to impedance data

Noise_Level=[0 0.001, 0.005,0.01]; % Define noise levels (0 = noise-free)
NL=Noise_Level(3); % Select a noise level; NL = 0 gives noise-free data
rng('default'); % Set random seed for reproducibility

% Add Gaussian noise proportional to the magnitude of Z_0
ReZ = real(Z_0)+ (NL*abs(Z_0).*randn(length(Z_0),1) + 0);
ImZ = imag(Z_0)+ (NL*abs(Z_0).*randn(length(Z_0),1) + 0);
Z_n=ReZ+1i.*ImZ; 

% Store frequency and impedance data
F = f(:); Z = Z_n(:);

%% Construction of the Loewner framework matrices and vectors 
% NOTE: Option `'real is set if one deals with data sets with imaginary
% entries
[L,Ls,V,W] = Loewner_Framework(1i*2*pi*F,Z,'real');

%% Determining the model order for SVD truncation 
% Estimate model order using rank of [L Ls] matrix (Loewner pencil)
k_LF=rank([L Ls]);

%% Construction of the reduced state-space model 
[Ek, Ak, Bk, Ck, s_LLs1] = state_space_mod(L,Ls,V,W,k_LF);

%% DRT calculation 
% Compute Distribution of Relaxation Times (DRT) from the model
[R_i_LF,tau_i_LF,pol_i_LF,res_i_LF] = DRT_values(Ek,Ak,Bk,Ck);

r=k_LF; % all extracted DRT elements are used in subsequent steps

%% Reproduction of impedance data using LF-derived DRT

% Intrapolation
Hr_DRT=[];
    for i1=1:length(F)
        hr=0;
    for j1=1:r
    a= R_i_LF(j1)./(1+1i*2*pi*F(i1)*tau_i_LF(j1));
        hr=hr+a;
    end
        Hr_DRT(i1,1)=hr;
    end

%Extrapolation
F_BIG = logspace(-2,6,8001); 
Hr_DRT_F_BIG=[];
    for i1=1:length(F_BIG)
        hr=0;
    for j1=1:r
    a= R_i_LF(j1)./(1+1i*2*pi*F_BIG(i1)*tau_i_LF(j1));
        hr=hr+a;
    end
        Hr_DRT_F_BIG(i1,1)=hr;
    end

%% Relative residuals

% Relative residuals between noisy data and reconstructed data 
ReZ_rel_DRT=((real(Z)-real(Hr_DRT))./(abs(Z)));
ImZ_rel_DRT=((imag(Z)-imag(Hr_DRT))./(abs(Z)));

% Relative difference between noise-free and noisy data
ReZ_rel_data=((real(Z_0)-real(Z))./(abs(Z_0)));
ImZ_rel_data=((imag(Z_0)-imag(Z))./(abs(Z_0)));

%% Plotting Figures

figure(1)
% Figure 1A: Nyquist plot
subplot('Position',[0.065 0.55 0.4 0.4])
plot(real(Z_n),-imag(Z_n),'o','color',[0 0 0],'LineWidth',3,...
    'markersize',7,'DisplayName','Noisy Data')
hold on
plot(real(Z_0),-imag(Z_0),'-','color',[0 1 1],'LineWidth',4,...
    'markersize',7,'DisplayName','Exact Data')
plot(real(Hr_DRT_F_BIG),-imag(Hr_DRT_F_BIG),'-','color',[1 0 0],...
    'markersize',7,'LineWidth',3,'DisplayName','LF')
grid on; axis equal;
xlabel('Re(Z) [\Omega]','FontSize',18,'FontWeight','normal')
ylabel('-Im(Z) [\Omega]','FontSize',18,'FontWeight','normal')
xlim([-20 220]); ylim([-5 95]);
xticks([0 40 80 120 160 200]);
xticklabels({'0','40','80','120','160','200'});
yticks([0 20 40 60 80]);
yticklabels({'0','20','40','60','80'});
legend('Location','NorthWest','FontSize',14,'FontWeight','normal',...
    'NumColumns',2)
annotation('textbox','String',{'A'},'FontSize',18,'FontWeight','bold',...
    'LineStyle','none','Position',[0.44,0.87,0.025,0.06]);

% Figure 1B: DRT plot
subplot('Position',[0.57 0.575 0.4 0.35])
stem(abs(tau_i_LF),abs(R_i_LF),'-.^','filled','color',[1 0 0],...
    'markersize',6,'LineWidth',2,'DisplayName','LF') 
set(gca,'xscal','log','yscal','log')
xlabel('\tau [s]','FontSize',18,'FontWeight','normal')
ylabel('|R_i| [\Omega]','FontSize',18,'FontWeight','normal')
grid on;
xlim([1e-13 1e3]); ylim([1e-4 1e4]);
xticks([1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1e0 1e2]);
xticklabels({'10^{-12}','10^{-10}','10^{-8}','10^{-6}','10^{-4}',...
    '10^{-2}','10^{0}','10^{2}'});
yticks([1e-4 1e-2 1e0 1e2 1e4]);
yticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}','10^{4}'});
legend('Location','NorthWest','FontSize',14,'FontWeight','normal',...
    'NumColumns',3)
annotation('textbox','String',{'B'},'FontSize',18,'FontWeight','bold',...
    'LineStyle','none','Position',[0.945,0.87,0.025,0.06]);

figure(2)
% Figure 2: Relative residuals
subplot('Position',[0.09,0.45,0.5,0.5])
plot(F,abs(ReZ_rel_data),':+','color',[0 0 0],...
    'markersize',7,'LineWidth',2,'DisplayName','Real-ND') 
hold on
plot(F,abs(ImZ_rel_data),'--o','color',[0 0 0],...
    'markersize',7,'LineWidth',2,'DisplayName','Imag-ND')
plot(F,abs(ReZ_rel_DRT),':+','color',[1 0 0],...
    'markersize',7,'LineWidth',2,'DisplayName','Real-LF') 
plot(F,abs(ImZ_rel_DRT),'--o','color',[1 0 0],...
    'markersize',7,'LineWidth',2,'DisplayName','Imag-LF')
set(gca,'xscal','log','yscal','log')
xlabel('Frequency [Hz]','FontSize',18,'FontWeight','normal')
ylabel('Residual [-]','FontSize',18,'FontWeight','normal')
grid on;
xlim([1e-2 1e6]); ylim([1e-17 1e5]);
xticks([1e-2 1e0 1e2 1e4 1e6]);
xticklabels({'10^{-2}','10^{0}','10^{2}','10^{4}','10^{6}'});
yticks([1e-15 1e-10 1e-5 1e0 1e5]);
yticklabels({'10^{-15}','10^{-10}','10^{-5}','10^{0}','10^{5}'});
legend('Location','NorthWest','FontSize',14,'FontWeight','normal',...
    'NumColumns',3)

figure(3)
% Figure 3: Singular value decay of the Loewner pencil ([L, Ls])
subplot('Position',[0.1 0.63 0.35 0.35])
semilogy(s_LLs1./s_LLs1(1),':ok','MarkerSize',7,'linewidth',2,...
    'DisplayName','SVs of [L, Ls]')
hold on
semilogy(k_LF,s_LLs1(k_LF)/s_LLs1(1),'^r','MarkerSize',6,'linewidth',...
    4,'DisplayName','Rank of [L, Ls]')
axis tight;grid on;
xlabel('Order of the model (r)')
ylabel('Singular value (SV)')
legend('Location','NorthEast','FontSize',14,'FontWeight',...
    'normal','NumColumns',1)

