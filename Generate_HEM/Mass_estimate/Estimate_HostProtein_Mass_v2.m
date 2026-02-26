%% Estimate_HostProtein_Mass_v2.m
% SynTwin  Host WT calibration (Generate_HEM / Mass_estimate)
%
% PURPOSE
%   Estimate an empirical relationship for host cell protein mass as a
%   function of the specific growth rate (mu), using literature data from
%   Bremer & Dennis (2008).
%
% CONTEXT IN SYNTWIN
%   This script is part of the internal calibration pipeline used to
%   reproduce the Wild-Type (WT) host parameterization reported in the
%   SynTwin publication.
%
%   The fitted host protein mass(mu) relationship contributes to the
%   baseline WT host description, which is later extended to the
%   Host Equivalent Model (HEM) under exogenous burden.
%
% DATA SOURCES (ATTRIBUTION)
%   Experimental values are taken from:
%     Bremer, H., & Dennis, P. P. (2008).
%     Modulation of chemical composition and other parameters of the cell
%     at different exponential growth rates.
%     EcoSal Plus, 3, 10.1128/ecosal.5.2.3.
%
%   Numerical values are included directly in this script as vectors for
%   reproducibility. See the repository-level
%   "Data Availability & Attribution Statement.md" for full attribution.
%
% INPUTS
%   None (literature vectors defined internally).
%
% OUTPUTS
%   - Fitted parameters for host protein mass(mu)
%   - Figures comparing model predictions to literature data
%
% USAGE
%   Run after initializing SynTwin paths:
%       ROOT = init_SynTwin();
%       run(fullfile(ROOT,'Generate_HEM','Mass_estimate',...
%           'Estimate_HostProtein_Mass_v2.m'));
%
% NOTES
%   - Growth rate mu units must be consistent ( 1/min;
%     conversions are handled explicitly inside the script).
%   - This script is included for reproducibility of paper results.
%     Typical SynTwin users do not need to run or modify it.
%
%
% -------------------------------------------------------------------------

% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin('bads',true);
% Ensure this script folder is on the path (for local functions in this folder)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);


td=[100 60, 40, 30, 24 20]; %(mins)
Mu_vector =log(2)./td';
% Experimental data:
% Protein dry mass (fg)
protein_Mass_old = [100 156 234 340 450 NaN]';
protein_Mass_new = [131 208 287 380 413 410]';
Exp_mh_vector = protein_Mass_new;
% Cell dry mass (fg)
cell_Mass_old = [148 258 433 641 865  NaN]';
cell_Mass_new = [226 374 555 774 921 1023]';    
Exp_cDW_vector = cell_Mass_new;
    
figure(1)
hold off
tiledlayout(1,3);
nexttile
plot(Mu_vector, protein_Mass_old, Mu_vector,protein_Mass_new,'Marker','*','MarkerSize',15,'LineWidth',2)
lgd = legend('$m_{h,\mathrm{old}}$','$m_{h,\mathrm{new}}$','Interpreter','latex','Location','northeast');
lgd.Layout.Tile = 'east';
lgd.FontSize = 16;
lgd.FontWeight = 'bold';
ylabel('$m_h$(fg)','Interpreter','latex'), xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex')
grid on
nexttile
plot(Mu_vector, cell_Mass_old, Mu_vector,cell_Mass_new,'Marker','*','MarkerSize',15,'LineWidth',2)
lgd = legend('$m_{cDW,\mathrm{old}}$','$m_{cDW,\mathrm{new}}$','Interpreter','latex','Location','northeast');
lgd.Layout.Tile = 'east';
lgd.FontSize = 16;
lgd.FontWeight = 'bold';
ylabel('$m_{cDW}$(fg)','Interpreter','latex'), xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex')
grid on
nexttile
plot(Mu_vector, cell_Mass_old./protein_Mass_old*100, Mu_vector,cell_Mass_new./protein_Mass_new*100,'Marker','*','MarkerSize',15,'LineWidth',2)
lgd = legend('$m_{cDW,old}/m_{h,old}$','$m_{cDW,new}/m_{h,new}$','Interpreter','latex','Location','northeast');
lgd.Layout.Tile = 'east';
lgd.FontSize = 16;
lgd.FontWeight = 'bold';
ylabel('$m_{cDW}/m_h$','Interpreter','latex'), xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex')
grid on

%% ESTIMATION OF HOST PROTEIN CONTENT MASS  mh

% We estimate a model of the form:
% MODEL  (Hill)
% predicted_mh = (p1+ 1e6*p2.*mu^p4)/(1+ 1000*p3.*mu^p4) 

options = bads('defaults');
options.Display='final';

x0 = [115,455,1030,3.3];   % Starting point for  p1, p2, p3, p4
lb = [100, 200, 500,2.5];  %lower expected bound for p1, p2, p3, p4
ub = [150,600,1500,4.0];  %upper expected bound for p1, p2, p3, p4 
plb = [110, 260, 500,3.0 ];  %Plausible lower bounds 
pub = [120,600,1500,3.5];  %Plausible upper bounds  

    % Run BADS, which returns the minimum X and its value FVAL.
    num_runs=10;   
    Results_mh=[];
parfor num_run=1:num_runs 
    J_mh=@(parameters) F_cost_mh_v2(parameters,Mu_vector,Exp_mh_vector);
    [params_mh, Jmin_value_mh] = bads(J_mh,x0,lb,ub,plb,pub,[],options)
    Results_mh=[Results_mh;[params_mh, Jmin_value_mh]]; 
end
format shorteng
params_found_mh = mean(Results_mh(:,1:4),1)
params_found_mh_std = std(Results_mh(:,1:4),0,1)

mu_test=0.002:0.0025:0.04;
% MODEL:
predicted_mh = ( params_found_mh(1) + 1e6*params_found_mh(2).*mu_test.^params_found_mh(4) )./( 1 +   1000*params_found_mh(3).*mu_test.^params_found_mh(4) );

figure(2)
hold off
tiledlayout(1,1);
nexttile
plot(Mu_vector,Exp_mh_vector,'Marker','square','MarkerSize',15)
hold on
plot(mu_test,predicted_mh,'Linewidth',3)
lgd = legend('$m_h$','$\hat{m}_h$','Interpreter','latex','Location','northeast');
lgd.Layout.Tile = 'east';
lgd.FontSize = 16;
lgd.FontWeight = 'bold';
ylabel('$m_h\, (\mathrm{fg})$','Interpreter','latex','Fontsize',16), xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex','Fontsize',16)
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
grid on

%%  ESTIMATION OF HOST CELL DRY WEIGHT  cDW

% We estimate a model of the form:
% MODEL  (Hill)
% predicted_cDW = (p1+ 1e6*p2.*mu^p4)/(1+ 1000*p3.*mu^p4) 

options = bads('defaults');
options.Display='final';

x0 = [2*115,2*455,1030,3.3];   % Starting point for  p1, p2, p3, p4
lb = [100, 200, 50,2.5];  %lower expected bound for p1, p2, p3, p4
ub = [2*150,2*600,2*1500,4.0];  %upper expected bound for p1, p2, p3, p4 
plb = [110, 260, 100,3.0 ];  %Plausible lower bounds 
pub = [2*120,2*600,2*1500,3.5];  %Plausible upper bounds  

    % Run BADS, which returns the minimum X and its value FVAL.
    num_runs=10;   
    Results_cDW=[];
parfor num_run=1:num_runs 
    J_cDW=@(parameters) F_cost_cDW_v2(parameters,Mu_vector,Exp_cDW_vector);
    [params_cDW, Jmin_value_cDW] = bads(J_cDW,x0,lb,ub,plb,pub,[],options)
    Results_cDW=[Results_cDW;[params_cDW, Jmin_value_cDW]]; 
end
format shorteng
params_found_cDW = mean(Results_cDW(:,1:4),1)
params_found_cDW_std = std(Results_cDW(:,1:4),0,1)

mu_test=0.002:0.0025:0.04;
% MODEL:
predicted_cDW = ( params_found_cDW(1) + 1e6*params_found_cDW(2).*mu_test.^params_found_cDW(4) )./( 1 +   1000*params_found_cDW(3).*mu_test.^params_found_cDW(4) );

figure(3)
hold off
tiledlayout(1,1);
nexttile
plot(Mu_vector,Exp_cDW_vector,'Marker','square','MarkerSize',15)
hold on
plot(mu_test,predicted_cDW,'Linewidth',3)
lgd = legend('$m_\mathrm{cDW}$','$\hat{m}_\mathrm{cDW}$','Interpreter','latex','Location','northeast');
lgd.Layout.Tile = 'east';
lgd.FontSize = 16;
lgd.FontWeight = 'bold';
ylabel('$m_\mathrm{cDW}\, (\mathrm{fg})$','Interpreter','latex','Fontsize',16), xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex','Fontsize',16)
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
grid on

%% PLOTING GLOBAL RESULTS
figure(3)
hold off
tiledlayout(1,3);
nexttile
plot(Mu_vector,Exp_mh_vector,'Marker','square','MarkerSize',10,'MarkerFaceColor','b')
hold on
plot(mu_test,predicted_mh,'Linewidth',3)
lgd = legend('$m_h$','$\hat{m}_h$','Interpreter','latex','Location','northeast');
lgd.Layout.Tile = 'east';
lgd.FontSize = 16;
lgd.FontWeight = 'bold';
ylabel('$m_h\, (\mathrm{fg})$','Interpreter','latex','Fontsize',16), xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex','Fontsize',16)
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
grid on

nexttile
plot(Mu_vector,Exp_cDW_vector,'Marker','square','MarkerSize',10,'MarkerFaceColor','b')
hold on
plot(mu_test,predicted_cDW,'Linewidth',3)
lgd = legend('$m_\mathrm{cDW}$','$\hat{m}_\mathrm{cDW}$','Interpreter','latex','Location','northeast');
lgd.Layout.Tile = 'east';
lgd.FontSize = 16;
lgd.FontWeight = 'bold';
ylabel('$m_\mathrm{cDW}\, (\mathrm{fg})$','Interpreter','latex','Fontsize',16), xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex','Fontsize',16)
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
grid on

nexttile
plot(Mu_vector,Exp_mh_vector./Exp_cDW_vector,'Marker','square','MarkerSize',10,'MarkerFaceColor','b')
hold on
plot(mu_test,predicted_mh./predicted_cDW,'Linewidth',3)
lgd = legend('$m_h/m_\mathrm{cDW}$','$\hat{m}_h/\hat{m}_\mathrm{cDW}$','Interpreter','latex','Location','northeast');
lgd.Layout.Tile = 'east';
lgd.FontSize = 16;
lgd.FontWeight = 'bold';
ylabel('$\frac{m_h}{m_\mathrm{cDW}}$','Interpreter','latex','Fontsize',16), xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex','Fontsize',16)
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
grid on

