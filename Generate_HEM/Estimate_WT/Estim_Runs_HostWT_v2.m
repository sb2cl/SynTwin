%% Estim_Runs_HostWT_v2
% SynTwin workflow script: parameter estimation of the *E. coli* wild-type (WT)
% host model at steady state using MEIGO (eSS).
%
% CONTEXT
%   This script performs the WT-host calibration that underpins the Host
%   Equivalent Model (HEM) used by SynTwin's host-aware digital twin.
%   The WT fit does not include exogenous burden parameters; it calibrates the
%   baseline host allocation model using literature data.
%
% DESCRIPTION
%   - Loads collated WT-host datasets (prepared as MATLAB structures) from
%     Bremer & Dennis (2008) and related sources.
%   - Runs global optimization with MEIGO's Enhanced Scatter Search (eSS)
%     to estimate steady-state host parameters.
%   - Saves the best parameter vectors and diagnostics for downstream use
%     (e.g., HEM surrogate generation).
%
% INPUTS
%   None (configuration is defined inside the script).
%
% OUTPUTS (saved to disk)
%   Results_Host_WT_v1.mat   (saved next to this script)
%
% DEPENDENCIES
%   - SynTwin_root.m and SynTwin_path.m for portable paths
%   - MATLAB/MEIGO/ (third-party optimizer)
%   - CostF_Host_WT_ss_v2.m (objective function)
%   - input_syn_data.m (loads/structures WT datasets)
%
% USAGE
%   Run from any folder:
%       Estim_Runs_HostWT_v2
%
% NOTES
%   This script is provided for reproducibility of the results reported in the
%   associated publication. Typical SynTwin users do not need to rerun it.
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
addpath(genpath(fullfile(ROOT,'Generate_HEM','Nut_estimate')));

% Global variables:
global model_p;

%% MODEL GENERAL PARAMETERS
model_p.mu_test = 0.005:0.0025:0.035;
model_p.mumax = 0.035;
model_p.f_s = model_p.mu_test/model_p.mumax;
model_p.nu_max = 20.5156*60; %min^{-1}

model_p.le_r = 24; 	%Ribosome occupancy length (aa)
model_p.lp_r = 195; %'Mean length of ribosomal proteins (aa)     
model_p.le_nr = 25; 	%Ribosome occupancy length (aa)
model_p.lp_nr = 333; %'Mean length of non-ribosomal proteins (aa)     
model_p.dm_r  =  0.16; %'Mean degradation rate of ribosomal mRNA (1/min)
model_p.dm_nr  =  0.2; %'Mean degradation rate of non-ribosomal mRNA (1/min)
model_p.N_r = 57;%56;% 57; %Number of protein types building up a ribosome (molec)			
model_p.N_nr = 1735; %Number of non ribosomal protein types being expressed at one given time (molec)
model_p.Em_r = model_p.lp_r/model_p.le_r*(1- (model_p.lp_r/(model_p.lp_r+model_p.le_r))^(model_p.lp_r/model_p.le_r)) ;  % Average Emr for ribosomal protein-coding genes
model_p.Em_nr =  model_p.lp_nr/model_p.le_nr*(1- (model_p.lp_nr/(model_p.lp_nr+model_p.le_nr))^(model_p.lp_nr/model_p.le_nr)) ;  % Average we
model_p.WEm_r = 1 + 1/model_p.Em_r;  % Average weight (1+1/Emr) for ribosomal protein-coding genes
model_p.WEm_nr =  1 + 1/model_p.Em_nr;  % Average weight (1+1/Emnr) for non-ribosomal protein-coding genes
model_p.m_aa = 182.6e-24; % Average amino acid mass (g/aa) 110 Da, 1Da = 1.6605e-24 g  
model_p.ribosome_mass = model_p.m_aa*model_p.N_r*model_p.lp_r; % 1e6 Da the protein content according to Bionumbers, which is ok with this estimation!
%Protein mass model m = ( 115.355e-15 +  455.799e-9*μ.^3.3448)./( 1 + 1028.37e3.*μ^3.3448 ) 
model_p.c0 = 115.355e-15;
model_p.c1 = 455.799e-9;
model_p.c2 = 1028.37e3;
model_p.c3 = 3.3448;
% nu_t_fs = @(f_s) ( model_p.nu_max*gamma1_fs*f_s./(gamma2_fs + f_s ) ); %aa/min
 model_p.gamma1_mu = 1.2629; 
 model_p.gamma2_mu =  0.0092;
 model_p.gamma2_fs = 0.2629;
 model_p.gamma1_fs = 1.2629;

 %% INITIAL ESTIMATIONS
model_p.ku_r =  132.68; %Dissotiation rate RBS-ribosome for ribosomal mRNA (1/min)
model_p.ku_nr =  6; %Dissotiation rate RBS-ribosome for non-ribosomal mRNA (1/min)            
model_p.kb_r =  7.82; %Assotiation rate RBS-ribosome for ribosomal mRNA (1/min/molec)
model_p.kb_nr =  15; %Assotiation rate RBS-ribosome for non-ribosomal mRNA (1/min/molec)                	
model_p.Omega_r = 2.4; % Average transcription rate for ribosomal proteins  (1/min)           
model_p.Omega_nr = 0.05; % Average transcription rate for non-ribosomal proteins  (1/min)        
model_p.k0_r =model_p.kb_r/model_p.ku_r;
model_p.k0_nr =model_p.kb_nr/model_p.ku_nr;
model_p.sigma0_r =model_p.nu_max/model_p.le_r/model_p.ku_r;
model_p.sigma0_nr =model_p.nu_max/model_p.le_nr/model_p.ku_nr;
model_p.Phi_m =  0.85; 
model_p.mu_profile = model_p.mu_test; %Initial guess of estimated mu profile
model_p.varphi_profile = 300*model_p.mumax*model_p.f_s; % %Initial guess of estimated varphi profile

 %% RUN OPTIMIZATIONS

num_runs = 50; 

Best_Results_Nr57=[];
for k=1:num_runs
     problem.f='CostF_Host_WT_ss_v2';    %script with the cost function to optimize
    % problem.x_L=[ 0.03,  2.0,      0.1,     20,     3.0,   0.04,     0.75]; % minimum expected values for k0_r, k0_nr, sigma0_r, sigma0_nr, omega_r, omega_nr, Phi_m
    % problem.x_U=[0.1,    20,           0.5,       80,     10,   0.14,   0.90]; %  maximum expected values
     problem.x_L=[0.06,  4,      0.10,       35,       4.05,   0.09,   0.8]; % minimum expected values for k0_r, k0_nr, sigma0_r, sigma0_nr, omega_r, omega_nr, Phi_m
    problem.x_U=[0.06,  10,      0.15,       70,       4.05,   0.09,   0.8]; %  maximum expected values

    opts.maxeval=1000; 
    %opts.maxtime=1e7;  
    opts.local.solver='fmincon';
    opts.ndiverse=90;
    Results=MEIGO(problem,opts,'ESS'); 
    vpa(Results.xbest')   
     Best_Results_Nr57= [Best_Results_Nr57, [Results.fbest, Results.xbest]'];
     k
end



%% STATISTICS OF ESTIMATED PARAMETERS

[const_index_sorted,cost_indexl_sorted_index]=sort(Best_Results_Nr57(1,:),'descend');
Best_Results_Nr57_sorted = Best_Results_Nr57(:,cost_indexl_sorted_index);
Best_Results_Nr57_sorted(1,:)= 1./Best_Results_Nr57_sorted(1,:);
% Best result is the last one.
% Get rid of  worst results. 
num_bad=round(0.6*size(Best_Results_Nr57_sorted,2));
Best_Results_Nr57_sorted_choice=Best_Results_Nr57_sorted(:,num_bad+1:end);
Weighted_sum_results_Nr57 = zeros(size(Best_Results_Nr57_sorted_choice,1),1);
Sum_weighs_results_Nr57 = 0;
for k=1:size(Best_Results_Nr57_sorted_choice,2)
    Weighted_sum_results_Nr57 = Weighted_sum_results_Nr57 + Best_Results_Nr57_sorted_choice(:,k)*Best_Results_Nr57_sorted_choice(1,k);
    Sum_weighs_results_Nr57 = Sum_weighs_results_Nr57 + Best_Results_Nr57_sorted_choice(1,k);
end
Weighted_average_results_Nr57= Weighted_sum_results_Nr57/Sum_weighs_results_Nr57
Regular_average_results_sorted_choice_Nr57=mean(Best_Results_Nr57_sorted_choice,2)
Regular_std_results_sorted_choice_Nr57=std(Best_Results_Nr57_sorted_choice,1,2)


%% EVALUATION OF THE AVERAGE ESTIMATION:

    estimated_m_params = Regular_average_results_sorted_choice_Nr57(2:8,1)';
    model_p.k0_r = estimated_m_params(1,1);
    model_p.k0_nr = estimated_m_params(1,2);
    model_p.sigma0_r = estimated_m_params(1,3);
    model_p.sigma0_nr = estimated_m_params(1,4);
    model_p.Omega_r = estimated_m_params(1,5);
    model_p.Omega_nr = estimated_m_params(1,6);
    model_p.Phi_m = estimated_m_params(1,7);

    f_s = model_p.f_s;

Eval_Host_results = eval_Host_WTss_fs_v2(model_p,f_s); 

%Protein mass model m = ( 115.355e-15 +  455.799e-9*μ.^3.3448)./( 1 + 1028.37e3.*μ^3.3448 ) 
mh = @(mu) ( model_p.c0 +  model_p.c1*mu.^model_p.c3)./( 1 + model_p.c2*mu.^model_p.c3 );    
% nu_t model 
 nu_t_fs = @(f_s) ( model_p.nu_max*model_p.gamma1_fs*f_s./(model_p.gamma2_fs + f_s ) ); %aa/min
% RBS strengths
KRBS_r = @(f_s) model_p.k0_r./(1 + model_p.sigma0_r/model_p.nu_max*nu_t_fs(f_s) );     
KRBS_nr = @(f_s) model_p.k0_nr./(1 + model_p.sigma0_nr/model_p.nu_max*nu_t_fs(f_s) );   
 % Resources recruitment strengths
 J_r = @(f_s,varphy) model_p.Em_r*model_p.Omega_r/model_p.dm_r*KRBS_r(f_s)./(1 + KRBS_r(f_s)/model_p.dm_r.*varphy);  
 J_nr = @(f_s,varphy) model_p.Em_nr*model_p.Omega_nr/model_p.dm_nr*KRBS_nr(f_s)./(1 + KRBS_nr(f_s)/model_p.dm_nr.*varphy);   
    
mu_test = model_p.mumax*f_s;
varphi_test = 0:0.05:35;
mp_estimated =  mh(f_s*model_p.mumax);
test_Jr = NaN*ones(length(f_s),length( varphi_test)); %Preallocate space
test_Jnr = NaN*ones(length(f_s),length( varphi_test)); %Preallocate space
test_WSum_h = NaN*ones(length(f_s),length( varphi_test)); %Preallocate space
test_Phi_t = NaN*ones(length(f_s),length( varphi_test));%Preallocate space
test_factor = NaN*ones(length(f_s),length( varphi_test)); %Preallocate space
mu_lines = NaN*ones(length(f_s),length( varphi_test)); %Preallocate space
varphi_lines = NaN*ones(length(f_s),length( varphi_test)); %Preallocate space
nut_vector= nu_t_fs(f_s);
for i=1:length(f_s)
    test_Jr(i,:) = J_r(f_s(i),varphi_test);
    test_Jnr(i,:)  = J_nr(f_s(i),varphi_test );
    test_WSum_h(i,:) = model_p.WEm_r *model_p.N_r *test_Jr(i,:) + model_p.WEm_nr *model_p.N_nr*test_Jnr(i,:);
    test_Phi_t(i,:) = model_p.N_r *test_Jr(i,:)./(1+test_WSum_h(i,:)); 
    test_factor(i,:)  = (model_p.Phi_m*test_Phi_t(i,:)).^2./(model_p.N_r *test_Jr(i,:)  + model_p.N_nr*test_Jnr(i,:) );
    varphi_lines(i,:)  = model_p.m_aa/model_p.ribosome_mass^2*test_factor(i,:)*mp_estimated(i).*nu_t_fs(f_s(i));
    mu_lines(i,:)  = 1/(model_p.N_r*model_p.lp_r)*model_p.Phi_m*test_Phi_t(i,:).*nu_t_fs(f_s(i));
 end % for f_s
    
 %% PLOT RESULTS

color = jet(length(f_s)); %one color for each substrate
modulus_factor = 2;

figure();
 tiledlayout(1,3);
nexttile
    k=0;
    hold on
 for i=1:length(f_s)
       if mod(i,modulus_factor)==0
             k=k+ (1-mod(i,modulus_factor));
             plot(varphi_test,varphi_lines(i,:),'Color',[color(i,:)],'LineWidth',3,'HandleVisibility','on')
             label{k} = strcat ('f(s)=',num2str( f_s(i) ) );
       else
           plot(varphi_test,varphi_lines(i,:),'Color',[color(i,:)],'LineWidth',3,'HandleVisibility','off')
       end
 end
  plot(varphi_test,varphi_test,'Color','k','LineWidth',4,'HandleVisibility','off')
grid on, ylabel('$\hat{\varphi}$','FontSize',18,'Interpreter','latex');
xlabel('$\varphi$','FontSize',18,'Interpreter','latex');
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.YAxis.FontSize = 18;
ax.YLabel.FontSize = 18;
nexttile
hold on 
 k=0;
    hold on
 for i=1:length(f_s)
       if mod(i,modulus_factor)==0
             k=k+ (1-mod(i,modulus_factor));
             plot(varphi_test,mu_lines(i,:),'Color',[color(i,:)],'LineWidth',3,'HandleVisibility','on')
             plot( Eval_Host_results.varphi_estimated, Eval_Host_results.mu_estimated,'linestyle','none','Marker','o','MarkerSize', 8,'MarkerFaceColor','k','HandleVisibility','off')
             label{k} = strcat ('f(s)=',num2str( f_s(i) ) );
       else
           plot(varphi_test,mu_lines(i,:),'Color',[color(i,:)],'LineWidth',3,'HandleVisibility','off')
       end
 end
grid on, ylabel('$\hat{\mu}$','FontSize',18,'Interpreter','latex');
xlabel('$\varphi$','FontSize',18,'Interpreter','latex');
 grid on
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.YAxis.FontSize = 18;
ax.YLabel.FontSize = 18;
nexttile
hold on 
 k=0;
    hold on
 for i=1:length(f_s)
       if mod(i,modulus_factor)==0
             k=k+ (1-mod(i,modulus_factor));
             plot(mu_test(i), Eval_Host_results.mu_estimated(i),'linestyle','none','Marker','o','MarkerSize', 12,'MarkerFaceColor',[color(i,:)],'MarkerEdgeColor',[color(i,:)],'HandleVisibility','on')
             label{k} = strcat ('f(s)=',num2str( f_s(i) ) );
       else
             plot(mu_test(i), Eval_Host_results.mu_estimated(i),'linestyle','none','Marker','o','MarkerSize', 12,'MarkerFaceColor',[color(i,:)],'MarkerEdgeColor',[color(i,:)],'HandleVisibility','off')
       end
 end
 plot(mu_test,mu_test,'Color','k','LineWidth',2,'HandleVisibility','off')
grid on, ylabel('$\hat{\mu}$','FontSize',18,'Interpreter','latex');
xlabel('$\mu$','FontSize',18,'Interpreter','latex');
 grid on
    lgd = legend(categorical(label),'Location','eastoutside');
    lgd.Layout.Tile = 'east';
    lgd.FontSize = 16;
    lgd.FontWeight = 'bold';
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.YAxis.FontSize = 18;
ax.YLabel.FontSize = 18;

%% SAVE RESULTS

file_WT_results = fullfile(SCRIPT_DIR,'Results_Host_WT_v2.mat');
save(file_WT_results,"Best_Results_OK","Best_Results_OK_sorted","Best_Results_OK_sorted_choice",...
    "Regular_average_results_sorted_choice_OK","Regular_std_results_sorted_choice_OK")