%% Estimate_Host_nu_v3.m
% SynTwin — Host WT calibration (Generate_HEM / Nut_estimate)
%
% PURPOSE
%   Fit an empirical relationship for the host peptide elongation rate
%   nu_t as a function of the specific growth rate mu, using literature
%   datasets compiled by Chure & Cremer (Flux Parity) and upstream sources
%   (incl. Bremer & Dennis).
%
% CONTEXT IN SYNTWIN
%   This script is part of the internal workflow used to reproduce the Host WT
%   calibration reported in the associated SynTwin paper. It provides a fitted
%   nu_t(mu) curve used downstream by the Host Equivalent Model (HEM) / digital twin.
%   Typical users of SynTwin do not need to run or modify this script.
%
% DATA SOURCES (ATTRIBUTION)
%   The datasets used here are derived from the Chure & Cremer "flux_parity"
%   repository (see the repository root file:
%       "Data Availability & Attribution Statement.md").
%   SynTwin may ship local Excel copies of the original CSVs for reproducibility.
%
% USAGE
%   Run from anywhere after initializing SynTwin paths:
%       ROOT = init_SynTwin('MEIGO', true);
%       run(fullfile(ROOT,'Generate_HEM','Nut_estimate','Estimate_Host_nu_v3.m'));
%
% INPUTS / OUTPUTS
%   Inputs : none (loads literature data via Load_Data_Nut_PhiR()).
%   Outputs: figures comparing PhiR(mu) and nu_t(mu) against data; the selected
%            fitted parameters are kept in the workspace.
%
% DEPENDENCIES
%   - init_SynTwin.m (SynTwin root)
%   - Load_Data_Nut_PhiR.m (Generate_HEM/Nut_estimate)
%   - MEIGO (if used for fitting) on MATLAB path (handled by init_SynTwin)
%   - hex2rgb (if not built-in in your MATLAB version, ensure it is available)
%
% NOTES
%   Units: mu is converted from 1/h to 1/min before fitting; nu_t is reported in aa/s.
%
% Copyright
%   SynTwin project. See repository LICENSE / README for terms of use.

clearvars;
close all;
dbstop if error
warning off;

% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin('MEIGO',true);
% Ensure this script folder is on the path (for local functions in this folder)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);
% Ensure  submodules used by this analysis are visible
addpath(genpath(fullfile(ROOT,'Generate_HEM','Third_party_data')));


Japan_5=['#00012a';'#fcc560';'#8fb1be';'#d53302';'#a0cbad'];
for i=1:5
     color(i,:) = hex2rgb(Japan_5(i,:));
end

global Mu_exp_nut Nu_t_exp Mu_exp_PhiR PhiR_exp


% Load data from Cremer 2023:

[Mu_exp_nut, Nu_t_exp, Mu_exp_PhiR, PhiR_exp, color_Cremer_phiR] =Load_Data_Nut_PhiR();
Mu_exp_nut = Mu_exp_nut/60;% to express in min^{-1}
Mu_exp_PhiR = Mu_exp_PhiR/60;% to express in min^{-1}


%Protein mass model m = ( 115.355e-15 +  455.799e-9*μ.^3.3448)./( 1 + 1028.37e3.*μ^3.3448 ) 
% with μ given in min^{-1}
model_p.c0 = 115.355e-15;
model_p.c1 = 455.799e-9;
model_p.c2 = 1028.37e3;
model_p.c3 = 3.3448;
mh = @(mu) ( model_p.c0 +  model_p.c1*mu.^model_p.c3)./( 1 + model_p.c2*mu.^model_p.c3 );


model_p.le_r = 24; 	%Ribosome occupancy length (aa)
model_p.lp_r = 195; %'Mean length of ribosomal proteins (aa)     
model_p.N_r = 57;%56;% 57; %Number of protein types building up a ribosome (molec)			
model_p.m_aa = 182.6e-24; % Average amino acid mass (g/aa) 110 Da, 1Da = 1.6605e-24 g  
model_p.ribosome_mass = model_p.m_aa*model_p.N_r*model_p.lp_r; % 1e6 Da the protein content according to Bionumbers, which is ok with this estimation!

PhiR_exp_base = PhiR_exp;
for i=1:30

    PhiR_exp = PhiR_exp_base/(0.95-i*0.01); % to account for mature ribosomes estimation

%% ESTIMATION OF THE FRACTION OF RIBOSOMES AS A FUNCTION OF MU

% We estimate a model PhiR=f(mu) of the form:
%
%     predicted_PhiR = a1 + a2*Mu_vector (mu_vector in min^{-1}

problem.x_L=[0.03, 3]; % lower expected bound for a1, a2
problem.x_U=[0.3,10]; %upper expected bound for a1, a2
problem.f='F_cost_PhiR';    %cost function to optimize
opts.maxeval=1000; 
opts.maxtime=1e7;  
opts.local.solver='fmincon';
opts.local.iterprint=1; 
Results_nut=MEIGO(problem,opts,'ESS');
 params_found_PhiR =Results_nut.xbest;
 a1(i) = params_found_PhiR(1)
 a2(i) = params_found_PhiR(2)


%% ASSOCIATED PREDICTION OF NU_t 

p1(i) = model_p.N_r*model_p.lp_r/a2(i)/60; %to give aa/sec
p2(i) = a1(i)/a2(i);
predicted_nut_mu = @(mu)  p1(i).*mu./( p2(i) + mu ); %Prediction in aa/sec

sqr_prediction_error = ( predicted_nut_mu(Mu_exp_nut) -Nu_t_exp ).^2  ;
J(i) = sum(sqr_prediction_error)/length(Nu_t_exp);

end

%% SELECTION OF "BEST " SOLUTION
 
mu_test=0.001:0.001:0.035;


[val,idx]= min(J)
factor = (0.95-idx*0.01)
predicted_PhiR_mu = @(mu)  a1(idx) + a2(idx).*mu; %Prediction 
p1_best = model_p.N_r*model_p.lp_r/a2(idx)/60; %to give aa/sec
p2_best = a1(idx)/a2(idx);
predicted_nut_mu = @(mu)  p1_best.*mu./( p2_best + mu ); %Prediction in aa/sec
PhiR_exp = PhiR_exp_base/(0.95-idx*0.01);

figure();
 tiledlayout(1,2);
nexttile
hold on 
plot(Mu_exp_PhiR, PhiR_exp,'linestyle','none','Marker','o','MarkerSize', 8,'MarkerFaceColor','k','HandleVisibility','on')
plot( mu_test, predicted_PhiR_mu( mu_test),'linestyle','-','Linewidth',4,'Color','r','HandleVisibility','on')
grid on, ylabel('$\hat{\Phi}_R$','FontSize',18,'Interpreter','latex');
xlabel('$\mu\, (\mathrm{min}^{-1})$','FontSize',18,'Interpreter','latex');
legend('$\Phi_R$','$\hat{\Phi}_R$','Interpreter','latex','Location','southeast','FontSize',16)
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.YAxis.FontSize = 18;
ax.YLabel.FontSize = 18;

nexttile
hold on 
plot(Mu_exp_nut, Nu_t_exp,'linestyle','none','Marker','o','MarkerSize', 8,'MarkerFaceColor','g','HandleVisibility','on')
plot( mu_test, predicted_nut_mu( mu_test),'linestyle','-','Linewidth',4,'Color','r','HandleVisibility','on')
grid on, ylabel('$\hat{\nu}_t$','FontSize',18,'Interpreter','latex');
xlabel('$\mu\, (\mathrm{min}^{-1})$','FontSize',18,'Interpreter','latex');
legend('$\nu_t$','$\hat{\nu}_t$','Interpreter','latex','Location','southeast','FontSize',16)
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.YAxis.FontSize = 18;
ax.YLabel.FontSize = 18;

