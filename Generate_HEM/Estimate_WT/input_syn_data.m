%% input_syn_data
% Load and assemble literature datasets used for WT-host model calibration.
%
% Used by: Estim_Runs_HostWT_v2, CostF_Host_WT_ss_v2, Analysis_WT_v1
%
% This script loads collated measurements (stored in MATLAB format for
% reproducibility) and exposes them as variables used by the WT-host
% estimation and validation routines.
%
% DATA SOURCES (upstream)
%   - Bremer, H., and Dennis, P. P. (2008). EcoSal Plus. DOI: 10.1128/ecosal.5.2.3
%   - Additional collations as documented in the associated paper/SI.
%
% FILES (in this SynTwin release)
%   - Generate_HEM/Third_party_data/Bremer_syn_data.mat
%
% NOTES
%   Typical SynTwin users do not need to modify this file.


global model_p;
global Bremer_syn_data;

load Bremer_syn_data.mat

model_p.mu_test = 0.005:0.005:0.035;
model_p.mumax = 0.035;
model_p.f_s = model_p.mu_test/model_p.mumax;
model_p.nu_max = 19.46*60; %min^{-1}
model_p.nu_max_aa_sec = 19.46; %aa*sec^{-1}

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

% Initial estimation of model parameters to be obtained:
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

%Protein mass model m = ( 115.355e-15 +  455.799e-9*μ.^3.3448)./( 1 + 1028.37e3.*μ^3.3448 ) 
model_p.c0 = 115.355e-15;
model_p.c1 = 455.799e-9;
model_p.c2 = 1028.37e3;
model_p.c3 = 3.3448;

 % nu_t_fs = @(f_s) ( model_p.nu_max*gamma1_fs*f_s./(gamma2_fs + f_s ) ); %aa/min
 model_p.gamma1_mu = 1.245; % 1.114;
 model_p.gamma2_mu =  0.00867; %0.00399103;
 model_p.gamma2_fs = model_p.gamma2_mu/model_p.mumax;
 model_p.gamma1_fs = 1.245; %1.114;

