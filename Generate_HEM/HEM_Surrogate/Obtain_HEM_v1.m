%% Obtain_HEM_v1
% SynTwin workflow script: generation of the Host Equivalent Model (HEM)
% surrogate used by the host-aware digital twin.
%
% PURPOSE
%   This script reproduces the HEM fitting pipeline used in the associated
%   publication and saves the resulting surrogate model as:
%       Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
%
%   The surrogate is later loaded by SynTwin estimation and post-processing
%   workflows (e.g., Get_synthesis_predictions / Get_synthesis_predictions_lite)
%   to compute synthesis-rate predictions as a function of growth rate and
%   circuit load.
%
% IMPORTANT NOTE (for users)
%   Typical SynTwin users do NOT need to run or modify this script. The
%   distributed file HEM_Surrogate.mat is provided for reproducibility of
%   the results reported in the paper.
%
% INPUTS
%   None (configuration is set within the script and required data files
%   are loaded from within the SynTwin repository).
%
% OUTPUTS (saved to disk)
%   - HEM_Surrogate.mat  (MATLAB struct containing the fitted HEM surrogate)
%     saved in: Generate_HEM/HEM_Surrogate/
%
% DEPENDENCIES
%   - init_SynTwin.m (repository root): initializes SynTwin paths.
%   - Generate_HEM/Estimate_WT/* : host fitting utilities called internally.
%   - Third-party biological datasets shipped with SynTwin for reproducibility
%     (see: "Data Availability & Attribution Statement.md" in repository root).
%
% USAGE
%   Obtain_HEM_v1
%
% SEE ALSO
%   eval_Host_fs_theta_v2, Get_synthesis_predictions, Get_synthesis_predictions_lite

clearvars; close all;
dbstop if error
warning on

% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin();  % Adds core folders (Scripts_base, Generate_HEM/HEM_Surrogate, etc.)

% Ensure this script folder is on the path (local helpers in this folder)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);

% Add Estimate_WT utilities (used internally by the HEM fitting pipeline)
addpath(SynTwin_path('Generate_HEM','Estimate_WT'));


model_p.le_r = 24; 	%Ribosome occupancy length (aa)
model_p.lp_r = 195; %'Mean length of ribosomal proteins (aa)     
model_p.le_nr = 25; 	%Ribosome occupancy length (aa)
model_p.lp_nr = 333; %'Mean length of non-ribosomal proteins (aa)     
model_p.dm_r  =  0.16; %'Mean degradation rate of ribosomal mRNA (1/min)
model_p.dm_nr  =  0.2; %'Mean degradation rate of non-ribosomal mRNA (1/min)
model_p.N_r = 57; %Number of protein types building up a ribosome (molec)			
model_p.N_nr = 1735; %Number of non ribosomal protein types being expressed at one given time (molec)
model_p.Em_r = model_p.lp_r/model_p.le_r*(1- (model_p.lp_r/(model_p.lp_r+model_p.le_r))^(model_p.lp_r/model_p.le_r)) ;  % Average Emr for ribosomal protein-coding genes
model_p.Em_nr =  model_p.lp_nr/model_p.le_nr*(1- (model_p.lp_nr/(model_p.lp_nr+model_p.le_nr))^(model_p.lp_nr/model_p.le_nr)) ;  % Average we
model_p.WEm_r = 1 + 1/model_p.Em_r;  % Average weight (1+1/Emr) for ribosomal protein-coding genes
model_p.WEm_nr =  1 + 1/model_p.Em_nr;  % Average weight (1+1/Emnr) for non-ribosomal protein-coding genes
model_p.m_aa = 182.6e-24; % Average amino acid mass (g/aa) 110 Da
Da = 1.6605e-24; %  Dalton unit g/Da
model_p.ribosome_mass = model_p.m_aa*model_p.N_r*model_p.lp_r; % 1e6 Da the protein content according to Bionumbers, which is ok with this estimation!

%Model parameters for WT
Results_WT_mean_values =[
 4.15  4.15  4.02  4.01  
0.09 0.09 0.09  0.09  
0.059 0.05 0.06 0.04  
9.23 10.31 11.55  6.60  
 0.14 0.10 0.12 0.15  
 63.78 73.65 71.67  66.84  
 0.79 0.79 0.80  0.81 ];
resul_us  = [3,4];
  model_p.Omega_r = mean(Results_WT_mean_values(1,resul_us)); %      
 model_p.Omega_nr = mean(Results_WT_mean_values(2,resul_us)); %   
 model_p.k0_r =  mean(Results_WT_mean_values(3,resul_us)); %
 model_p.k0_nr = mean(Results_WT_mean_values(4,resul_us)); % 
 model_p.sigma0_r = mean(Results_WT_mean_values(5,resul_us)); % 
 model_p.sigma0_nr = mean(Results_WT_mean_values(6,resul_us)); % 
 model_p.Phi_m = mean(Results_WT_mean_values(7,resul_us)); % 

 %Protein mass model m = ( 115.355e-15 +  455.799e-9*μ.^3.3448)./( 1 + 1028.37e3.*μ^3.3448 ) 
model_p.c0 = 115.355e-15;
model_p.c1 = 455.799e-9;
model_p.c2 = 1028.37e3;
model_p.c3 = 3.3448;

%nu_t model nu_t = @(f_s) ( model_p.nu_max*gamma1*f_s./(gammap2 + f_s ) ); %aa/min
model_p.nu_max = 20.5156*60; %min^{-1}
model_p.gamma1_fs = 1.2629;
model_p.gamma2_fs = 0.2629;


 %%%%% Functions for the Host WT cell
    %Protein mass model m = ( 115.355e-15 +  455.799e-9*μ.^3.3448)./( 1 + 1028.37e3.*μ^3.3448 ) 
         mh = @(mu) ( model_p.c0 +  model_p.c1*mu.^model_p.c3)./( 1 + model_p.c2*mu.^model_p.c3 );
    % nu_t model as a function of normalised substrate f(s) = mu/mu_max
         nu_t = @(f_s) ( model_p.nu_max*model_p.gamma1_fs*f_s./(model_p.gamma2_fs + f_s ) ); %aa/min
    % RBS strengths
         KRBS_r = @(f_s) model_p.k0_r./(1 + model_p.sigma0_r/model_p.nu_max*nu_t(f_s) );     
         KRBS_nr = @(f_s) model_p.k0_nr./(1 + model_p.sigma0_nr/model_p.nu_max*nu_t(f_s) );   
     % Resources recruitment strengths
        J_r = @(f_s,varphy) model_p.Em_r*model_p.Omega_r/model_p.dm_r*KRBS_r(f_s)./(1 + KRBS_r(f_s)/model_p.dm_r.*varphy);  
        J_nr = @(f_s,varphy) model_p.Em_nr*model_p.Omega_nr/model_p.dm_nr*KRBS_nr(f_s)./(1 + KRBS_nr(f_s)/model_p.dm_nr.*varphy);     

 %%%% External substrate function as the function f(s) = s/(K_s+s) = mu/mu_max
        f_s = 0.005:0.005:1; %this contains the fraction of substrate f(s). This is the external input to the model
        model_p.f_s = f_s;

%%%% Initial estimation of the mu and varphi profiles (required for the
%%%% iterative procedure in eval_Hostss_fs_theta_v1.m
    model_p.mu_max = 0.0344;
    model_p.mu_test = model_p.mu_max*f_s;
    model_p.mu_profile =  model_p.mu_test; %Initial guess of estimated mu profile
    model_p.varphi_profile = 700*model_p.mu_max*f_s; % %Initial guess of estimated varphi profile

%%%% Exogenous loading (theta=0 corresponds to the wild-type)
   % Theta_values = [0,10,100,250, 500, 750, 1000];
    %Theta_values = 0:0.05:500; 
    %Theta_values = logspace(-2,3,750); 
    Theta_values = logspace(-4,3.5,800); 

%%%% NOW we simulate the steady state %%%%%%%%%%%
Matrix_varphi_estimated = [];
Matrix_mu_estimated = [];
Matrix_r_estimated = [];
Matrix_rT_estimated = [];
Matrix_rt_estimated = [];
Matrix_ra_estimated = [];
Matrix_estimated_Jr = [];
Matrix_estimated_Jnr = [];
for i=1:length(Theta_values)
    Eval_Host_theta_results= eval_Host_fs_theta_v2(model_p, f_s,Theta_values(i)); 
    Matrix_varphi_estimated = [Matrix_varphi_estimated; Eval_Host_theta_results.varphi_estimated];
    Matrix_mu_estimated = [Matrix_mu_estimated; Eval_Host_theta_results.mu_estimated];
    Matrix_estimated_Jr = [Matrix_estimated_Jr; Eval_Host_theta_results.J_r];
    Matrix_estimated_Jnr = [Matrix_estimated_Jnr; Eval_Host_theta_results.J_nr];
    Matrix_r_estimated = [Matrix_r_estimated; Eval_Host_theta_results.free_r];
    Matrix_rT_estimated = [Matrix_rT_estimated; Eval_Host_theta_results.total_r_T];
    Matrix_rt_estimated = [Matrix_rt_estimated; Eval_Host_theta_results.active_r_t];
    Matrix_ra_estimated = [Matrix_ra_estimated; Eval_Host_theta_results.available_r_a];
end

%%%% Evaluation of the results
   % Load data from Cremer 2023:
[Mu_exp_nut, Nu_t_exp, Mu_exp_PhiR, PhiR_exp, color_Cremer_phiR] = Load_Data_Nut_PhiR();
Mu_exp_nut = Mu_exp_nut/60;% to express in min^{-1}
Mu_exp_PhiR = Mu_exp_PhiR/60;% to express in min^{-1}

   
F_s_values = f_s; %for saving purposes

% Build a struct with the results for the Host Equivalent Model

HEM.F_s_values = F_s_values;
HEM.Theta_values = Theta_values;
HEM.Matrix_varphi_values = Matrix_varphi_estimated;
HEM.Matrix_mu_values = Matrix_mu_estimated;
HEM.Matrix_r_values = Matrix_r_estimated;

file_HEM = fullfile(SCRIPT_DIR,'HEM_Surrogate.mat');
save(file_HEM,'HEM','-mat');

figures_plot=true;
if figures_plot==true
     color = jet(length(Theta_values)); %one color for each exogenous loading
     modulus_factor = 50;
    f1 = figure(1);
    hold off
    tiledlayout(1,2);
    nexttile
%    label={}; 
     k=0;
      hold on
    for i=1:length(Theta_values)
         if mod(i,modulus_factor)==0
             k=k+ (1-mod(i,modulus_factor));
             plot(f_s,Matrix_mu_estimated(i,:),'Color',[color(i,:)],'LineWidth',2,'HandleVisibility','on')
       %      label{k} = strcat ('Load:',num2str( Theta_values(i) )  );
               else
           plot(f_s,Matrix_mu_estimated(i,:),'Color',[color(i,:)],'LineWidth',2,'HandleVisibility','off')
       end
    end
    ax = gca;
    ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
ax.XLabel.FontSize = 20;
ax.YLabel.FontSize = 20;
    hold off
    ylabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex'), xlabel('$f(s_n)$','Interpreter','latex')
    grid on
%    lgd = legend(categorical(label),'Location','eastoutside');
    lgd.Layout.Tile = 'east';
    lgd.FontSize = 18;
    lgd.FontWeight = 'bold';
nexttile
    label={}; 
     k=0;
      hold on
    for i=1:length(Theta_values)
         if mod(i,modulus_factor)==0
             k=k+ (1-mod(i,modulus_factor));
             plot(f_s,Matrix_varphi_estimated(i,:),'Color',[color(i,:)],'LineWidth',2,'HandleVisibility','on')
          %   label{k} = strcat ('Load:',num2str( Theta_values(i), '%3.2f' )  );
             label{k} = num2str( Theta_values(i), '%3.2f' );
               else
           plot(f_s,Matrix_varphi_estimated(i,:),'Color',[color(i,:)],'LineWidth',2,'HandleVisibility','off')
       end
    end
    ax = gca;
    ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
ax.XLabel.FontSize = 20;
ax.YLabel.FontSize = 20;
    hold off
    ylabel('$\varphi\, (\mathrm{molec}\cdot\mathrm{min}^{-1})$','Interpreter','latex'), xlabel('$f(s_n)$','Interpreter','latex')
    grid on
    lgd = legend(categorical(label),'Location','eastoutside');
    lgd.Layout.Tile = 'east';
    lgd.FontSize = 18;
    lgd.FontWeight = 'bold';
    title(lgd,'Load ($\theta$)','Interpreter','latex')
   

f2= figure(2);
tiledlayout(1,1);
nexttile
 label={}; 
 hold on
    for i=1:length(Theta_values)
         if mod(i,modulus_factor)==0
             k=k+ (1-mod(i,modulus_factor));
             plot(f_s,Matrix_varphi_estimated(i,:),'Color',[color(i,:)],'LineWidth',2,'HandleVisibility','on')
               label{k} = strcat ('Load:',num2str( Theta_values(i) )  );
               else
            plot(f_s,Matrix_varphi_estimated(i,:),'Color',[color(i,:)],'LineWidth',2,'HandleVisibility','off')
       end
    end
    ax = gca;
    ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;
    ylabel('$\varphi$','Interpreter','latex'), xlabel('$f(s)$','Interpreter','latex')
    grid on
    lgd = legend(categorical(label),'Location','eastoutside');
    lgd.Layout.Tile = 'east';
    lgd.FontSize = 16;
    lgd.FontWeight = 'bold';
    hold off

f3= figure(3);
tiledlayout(1,2);
nexttile
    for i=1:length(Theta_values)
        plot(Matrix_varphi_estimated(i,:),f_s,'Color',[color(i,:)],'LineWidth',2,'HandleVisibility','on')
        hold on
    end
    ax = gca;
    ax.FontSize = 14;
    xlabel('$\varphi$','Interpreter','latex'), ylabel('$f(s)$','Interpreter','latex')
    grid on
nexttile
 label={}; 
         hold on

    for i=1:length(Theta_values)
         if mod(i,modulus_factor)==0
             k=k+ (1-mod(i,modulus_factor));
              plot(Matrix_estimated_Jnr(i,:),f_s,'Color',[color(i,:)],'LineWidth',2,'HandleVisibility','on')
         label{k} = strcat ('Load:',num2str( Theta_values(i) )  );
         else
            plot(Matrix_estimated_Jnr(i,:),f_s,'Color',[color(i,:)],'LineWidth',2,'HandleVisibility','off')
       end
    end
    ax = gca;
    ax.FontSize = 14;
    xlabel('$J_{nr}$','Interpreter','latex'), ylabel('$f(s)$','Interpreter','latex')
    grid on
    lgd = legend(categorical(label),'Location','eastoutside');
    lgd.Layout.Tile = 'east';
    lgd.FontSize = 16;
    lgd.FontWeight = 'bold';
    hold off

f4= figure(4);
tiledlayout(1,1);
nexttile
 label={}; 
  hold on
    for i=1:length(Theta_values)
         if mod(i,modulus_factor)==0
             k=k+ (1-mod(i,modulus_factor));
              plot(Matrix_varphi_estimated(i,:),Matrix_mu_estimated(i,:),'Color',[color(i,:)],'LineWidth',2,'HandleVisibility','on')
               label{k} = strcat ('Load:',num2str( Theta_values(i) )  );
         else
             plot(Matrix_varphi_estimated(i,:),Matrix_mu_estimated(i,:),'Color',[color(i,:)],'LineWidth',2,'HandleVisibility','off')
         end
    end
    ax = gca;
    ax.FontSize = 14;
    xlabel('$\varphi$','Interpreter','latex'), ylabel('$f(s)$','Interpreter','latex')
    grid on
    lgd = legend(categorical(label),'Location','eastoutside');
    lgd.Layout.Tile = 'east';
    lgd.FontSize = 16;
    lgd.FontWeight = 'bold';
    hold off

end %figures plot

