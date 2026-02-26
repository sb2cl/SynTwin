%% Analysis_WT_v1
% SynTwin workflow script: validation and visualization of the fitted WT-host
% model against literature datasets.
%
% Used by: (standalone analysis script)
%
% This script loads previously estimated WT-host parameters (Results_Host_WT_v1)
% and compares model predictions (e.g., ribosomal mass fraction relationships)
% to collated measurements used for calibration/validation.
%
% OUTPUT
%   Generates figures for the software/paper supplementary material.
%
% NOTES
%   Provided for reproducibility of the results reported in the associated
%   publication. Typical SynTwin users do not need to rerun it.
%
% Copyright
%   SynTwin project. See repository LICENSE / README for terms of use.

clearvars
close all

% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin('MEIGO',true);
% Ensure this script folder is on the path (for local functions in this folder)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);
% Ensure  submodules used by this analysis are visible
addpath(genpath(fullfile(ROOT,'Generate_HEM','Third_party_data')));
addpath(genpath(fullfile(ROOT,'Generate_HEM','Nut_estimate')));

%% DECLARATION OF PARAMETERS AND FUNCTIONS

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

model_p.mumax = 0.035;
model_p.nu_max = 20.5156*60; %min^{-1}
model_p.nu_max_aa_sec = 20.5156; %aa*sec^{-1}

mu_test = 0.002:0.001:0.035;
f_s = mu_test/model_p.mumax;
model_p.mu_profile = mu_test; %Initial guess of estimated mu profile
model_p.varphi_profile = 300*model_p.mumax*f_s; % %Initial guess of estimated varphi profile

% The 5th result is very sensitive 
Results_WT_mean_values =[
 4.15  4.15  4.02  4.01 3.82  
0.09 0.09 0.09  0.09 0.09 
0.059 0.05 0.06 0.04 0.04 
9.23 10.31 11.55  6.60 7.76 
 0.14 0.10 0.12 0.15 0.13 
 63.78 73.65 71.67  66.84 66.78 
 0.79 0.79 0.80  0.81 0.85];
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
mh = @(mu) ( model_p.c0 +  model_p.c1*mu.^model_p.c3)./( 1 + model_p.c2*mu.^model_p.c3 );    

% nu_t model 
 model_p.gamma1_mu = 1.2629; 
 model_p.gamma2_mu =  0.0092;
 model_p.gamma2_fs = 0.2629;
 model_p.gamma1_fs = 1.2629;
nu_t_fs = @(f_s) ( model_p.nu_max*model_p.gamma1_fs*f_s./(model_p.gamma2_fs + f_s ) ); %aa/min
 

% RBS strengths
KRBS_r = @(f_s) model_p.k0_r./(1 + model_p.sigma0_r/model_p.nu_max*nu_t_fs(f_s) );     
KRBS_nr = @(f_s) model_p.k0_nr./(1 + model_p.sigma0_nr/model_p.nu_max*nu_t_fs(f_s) );   
 % Resources recruitment strengths
 J_r = @(f_s,varphy) model_p.Em_r*model_p.Omega_r/model_p.dm_r*KRBS_r(f_s)./(1 + KRBS_r(f_s)/model_p.dm_r.*varphy);  
 J_nr = @(f_s,varphy) model_p.Em_nr*model_p.Omega_nr/model_p.dm_nr*KRBS_nr(f_s)./(1 + KRBS_nr(f_s)/model_p.dm_nr.*varphy);   
    
%% GETTING MU AND VARPHI FOR FS and GENERIC MU&VARPHI "LINES"

Eval_Host_results = eval_Host_WTss_fs_v2(model_p,f_s); 

mu_test = model_p.mumax*f_s;
varphi_test = 0:0.05:30;
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

%% GETTING ESTIMATED RBS STRENGHS, RRSs, NUMBER OF RIBOSOMES AND FRACTIONS OF RIBOSOMES

 % WE HAVE:
% Eval_Host_results.varphi_estimated 
% Eval_Host_results.mu_estimated 
% Eval_Host_results.WSum_h 
% Eval_Host_results.Sum_h = N_r *Jr+ N_nr*Jnr;
% Eval_Host_results.Phi_r = fraction of ribosomes as Nr*Jr/Sum_h
% Eval_Host_results.J_r 
% Eval_Host_results.J_nr 
% Eval_Host_results.nu_t_fs 
% Eval_Host_results.fs 
% Eval_Host_results.KRBS_r 
% Eval_Host_results.KRBS_nr 

% NOTICE: 
% free ribosomes r=varphi/mu
free_r = Eval_Host_results.varphi_estimated./Eval_Host_results.mu_estimated;
% available mature ribosomes   r_a = r*(1+ WSum_h) = varphi/mu*(1+ WSum_h)
available_r_a = (1+ Eval_Host_results.Sum_h).*free_r ;
% bound ribosomes r_b = WSum_h/(1+ WSum_h)*r_a = WSum_h*varphi/mu =
% WSum_h*r
bound_r_b = Eval_Host_results.Sum_h.*free_r;
% active translating ribosomes r_t = Sum_h/(1+ WSum_h)*r_a = Sum_h*r
active_r_t = Eval_Host_results.Sum_h.*free_r;
estimated_Phi_r = model_p.ribosome_mass*active_r_t./mp_estimated;
%  total ribosomes r_T = r_a/Phi_m 
total_r_T = available_r_a/model_p.Phi_m;

%% LOAD DATA FROM CREMER'S LAB

[Mu_exp_nut, Nu_t_exp, Mu_exp_PhiR, PhiR_exp, color_Cremer_phiR] =Load_Data_Nut_PhiR();
Mu_exp_nut = Mu_exp_nut/60;% to express in min^{-1}
Mu_exp_PhiR = Mu_exp_PhiR/60;% to express in min^{-1}


%% PLOTS

color = jet(length(f_s)); %one color for each substrate
modulus_factor = 4;
%figure 1
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
ax.YLabel.FontSize = 18;
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
ax.XLabel.FontSize = 18;
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

%figure 2
figure();
 tiledlayout(1,2);
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

%figure 3

mu_test=0.001:0.001:0.035;
a1= 0.0659744;
a2 = 7.14928;
predicted_PhiR_mu = @(mu)  a1 + a2.*mu; %Prediction from nu_t

figure();
 tiledlayout(1,2);
nexttile
hold on 
plot(Mu_exp_PhiR, PhiR_exp/model_p.Phi_m,'linestyle','none','Marker','o','MarkerSize', 8,'MarkerFaceColor','k','HandleVisibility','on')
plot( mu_test, predicted_PhiR_mu( mu_test),'linestyle','-','Linewidth',2,'Color','g','HandleVisibility','off')
plot( Eval_Host_results.mu_estimated, estimated_Phi_r,'linestyle','-','Linewidth',4,'Color','r','HandleVisibility','on')
grid on, ylabel('$\Phi_R$','FontSize',18,'Interpreter','latex');
xlabel('$\mu\, (\mathrm{min}^{-1})$','FontSize',18,'Interpreter','latex');
legend('$\Phi_R$','$\hat{\Phi}_R$','Interpreter','latex','Location','southeast','FontSize',16)
 grid on
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.YAxis.FontSize = 18;
ax.YLabel.FontSize = 18;
nexttile
hold on 
R_s=[bound_r_b ;free_r;(total_r_T-available_r_a)];
area(Eval_Host_results.mu_estimated,R_s','LineWidth',0.1)
plot( Eval_Host_results.mu_estimated, total_r_T,'linestyle','-','Linewidth',4,'Color','k','HandleVisibility','on')
grid on, ylabel('$r_b,r,r_i$','FontSize',18,'Interpreter','latex');
xlabel('$\mu\, (\mathrm{min}^{-1})$','FontSize',18,'Interpreter','latex');
legend('Bound ribosomes', 'Free ribosomes','Immature ribosomes', 'Total ribosomes', 'Location','northwest','FontSize',16)
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.YAxis.FontSize = 18;
ax.YLabel.FontSize = 18;
