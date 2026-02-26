function Pi_pred_values = Get_synthesis_predictions_lite(HEM,model_c, Mu_vector,RBS_k0_sigma0,RBS_inv_sigma0,Omega,Gen_cn)
%% Get_synthesis_predictions_lite
% SynTwin core function: lightweight digital-twin prediction of TU synthesis rate.
%
% PURPOSE
%   Compute the predicted synthesis-rate trajectory Pi (Π) for an exogenous
%   Transcriptional Unit (TU) over a provided growth-rate vector Mu_vector,
%   using the Host Equivalent Model (HEM) digital twin.
%
%   This "lite" version is intentionally minimal and fast: it returns only
%   the predicted synthesis-rate values (Pi_pred_values) and avoids computing
%   auxiliary host variables and parametric sensitivities. It is therefore
%   suitable for computationally demanding optimization loops (e.g., BADS/eSS).
%
% MODEL CONTEXT
%   SynTwin predicts TU synthesis by coupling:
%     (i)  a host-aware HEM surrogate (resource allocation / physiology), and
%     (ii) a mechanistic TU synthesis model parameterized by:
%          - Gen_cn        : effective gene copy number
%          - Omega         : promoter strength (transcriptional capacity)
%          - RBS_k0_sigma0 : RBS intrinsic initiation capacity (IIC; k0/sigma0)
%          - RBS_inv_sigma0: RBS sensitivity-related parameter (1/sigma0)
%
%   The digital twin uses Mu_vector as an exogenous input to constrain the host
%   state consistently with the experimentally observed growth regime.
%
% INPUTS
%   HEM             : struct containing the HEM surrogate data (loaded from
%                     Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat).
%   model_c         : struct with circuit/CDS constants (e.g., lp_c, le_c, dm_c,
%                     Em_c, WEm_c, etc.).
%   Mu_vector       : (n x 1) vector of growth-rate values (1/min).
%   RBS_k0_sigma0   : scalar, intrinsic initiation capacity (k0/sigma0) of the RBS.
%   RBS_inv_sigma0  : scalar, inverse sensitivity parameter (1/sigma0).
%   Omega           : scalar, promoter strength parameter.
%   Gen_cn          : scalar, effective gene copy number.
%
% OUTPUTS
%   Pi_pred_values  : (n x 1) vector of predicted synthesis-rate values Π over Mu_vector.
%
% DEPENDENCIES
%   - SynTwin initialization is recommended before use:
%       ROOT = init_SynTwin(...);
%   - Requires the HEM surrogate struct (HEM) generated/provided with SynTwin.
%   - Shared mathematical subroutines are contained within SynTwin core folders
%     (added by init_SynTwin).
%
% USAGE
%   Pi_pred = Get_synthesis_predictions_lite(HEM, model_c, Mu_vec, ...
%              RBS_k0_sigma0, RBS_inv_sigma0, Omega, Gen_cn);
%
% NOTES
%   - For the full set of outputs (host variables, internal coefficients, and
%     relative sensitivities), use Get_synthesis_predictions.m.
%   - This function is designed for repeated calls inside optimizers; keep
%     Mu_vector reasonably sized for speed in large-scale estimation runs.
%
% AUTHORS / VERSION
%   SynTwin development team. This release corresponds to the implementation
%   used to generate the computational results reported in the associated paper.
 
 % nu_t model as a function of normalised substrate f(s) = mu/mu_max
    model_p.nu_max = 20.5156*60; %min^{-1}
    model_p.gamma1_fs = 1.2629;
    model_p.gamma2_fs = 0.2629;
    nu_t = @(f_s) ( model_p.nu_max*model_p.gamma1_fs*f_s./(model_p.gamma2_fs + f_s ) ); %aa/min
 % RBS strength
     KRBS_c= @(f_s, k0_sigma0, inv_sigma0) k0_sigma0./(inv_sigma0 + 1/model_p.nu_max*nu_t(f_s) );    
  % Resources recruitment strength
     J_c = @(f_s,varphy,k0_sigma0, inv_sigma0,Omega) model_c.Em_c*Omega/model_c.dm_c*KRBS_c(f_s, k0_sigma0, inv_sigma0)./(1 + KRBS_c(f_s, k0_sigma0, inv_sigma0)/model_c.dm_c.*varphy);  
  % Loading function
     Load_theta = @(f_s,varphy,k0_sigma0, inv_sigma0,Omega,Nc) model_c.WEm_c*Nc* J_c(f_s,varphy,k0_sigma0, inv_sigma0,Omega) ;
  % Synthesis rate
     Pi_syn = @(f_s,varphy,mu,k0_sigma0, inv_sigma0,Omega,Nc) Nc/model_c.lp_c*J_c(f_s,varphy,k0_sigma0, inv_sigma0,Omega).*varphy./mu.*nu_t(f_s);
 Pi_pred_values = zeros(length(Mu_vector),1);
 for r=1:length(Mu_vector)
         % For  mu_exp(r) we find the corresponding flux varphi in the HEM  for each potential value of the substrate function f_s: 
            for s=1:length(HEM.F_s_values)
                HEM_varphi{s}.values = HEM.Matrix_varphi_values(:,s); 
                HEM_mu{s}.values = HEM.Matrix_mu_values(:,s); 
                Current_TU.Varphi_exp_avg(s) = Find_varphi(HEM_mu{s}.values,HEM_varphi{s}.values, Mu_vector(r) );
           % For varphi  we find the corresponding load induced by the circuit:
                Current_TU.Load_exp_avg(s) = Load_theta(HEM.F_s_values(s),Current_TU.Varphi_exp_avg(s), RBS_k0_sigma0,RBS_inv_sigma0,Omega,Gen_cn);
             % For that load we find the flux in HEM:
                 Current_TU.Varphi_HEM_avg(s)  = Find_varphiHEM_thetaExo(HEM_varphi{s}.values,HEM.Theta_values, Current_TU.Load_exp_avg(s) ); 
           % For that flux varphi  we find the corresponding synthesis rate:
                Current_TU.Pi_pred_avg(s) =  Pi_syn(HEM.F_s_values(s),Current_TU.Varphi_exp_avg(s), Mu_vector(r), RBS_k0_sigma0,RBS_inv_sigma0,Omega,Gen_cn);
            end
        % The intersection point between Current_TU.Varphi_exp_avg and Current_TU.varphi_HEM for some value of the substrate 
        % function f_s,  gives the actual flux of resources
        % (a single reconciliated value) and the corresponding load at the single common value of f_s corresponding to that experimentyal mu_avg.
        [Current_TU.varphi_common_avg, Current_TU.theta_common_avg, Current_TU.fs_common_avg]  = Find_HEMLoad_varphi_fs_common(Current_TU.Varphi_exp_avg,...
                                                                                                                                                                                Current_TU.Varphi_HEM_avg, Current_TU.Load_exp_avg, HEM.F_s_values);
        Current_TU.Pi_pred_common_avg =  Pi_syn(Current_TU.fs_common_avg,Current_TU.varphi_common_avg, Mu_vector(r),RBS_k0_sigma0,RBS_inv_sigma0,Omega,Gen_cn);
        % We save the "reconciliated" results:
        Pi_pred_values(r,1) = Current_TU.Pi_pred_common_avg;
end  % for values in Mu_exp_vector

end