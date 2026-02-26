function Synthesis_predictions = Get_synthesis_predictions(HEM,model_c, Mu_vector,RBS_k0_sigma0,RBS_inv_sigma0,Omega,Gen_cn)
%% Get_synthesis_predictions
% SynTwin core function: digital-twin prediction of TU synthesis and sensitivities.
%
% PURPOSE
%   Compute host-aware digital-twin predictions for an exogenous Transcriptional
%   Unit (TU) over a provided growth-rate vector Mu_vector. In addition to the
%   predicted synthesis rate Π, this full version returns internal host/circuit
%   variables and relative parametric sensitivities needed for:
%     - Results-tensor construction,
%     - Sensitivity mapping at experimental operating points,
%     - Monte Carlo uncertainty propagation,
%     - Practical identifiability analysis (Jacobian-based diagnostics).
%
% MODEL CONTEXT
%   SynTwin predicts TU synthesis by coupling:
%     (i)  a Host Equivalent Model (HEM) surrogate capturing resource allocation
%          and physiology, and
%     (ii) a mechanistic TU synthesis model parameterized by biopart quantities:
%          - Gen_cn        : effective gene copy number
%          - Omega         : promoter strength (transcriptional capacity)
%          - RBS_k0_sigma0 : RBS intrinsic initiation capacity (IIC; k0/sigma0)
%          - RBS_inv_sigma0: RBS sensitivity-related parameter (1/sigma0)
%
%   The digital twin uses Mu_vector as an exogenous input that constrains the
%   host state, enabling context-aware predictions consistent with measured
%   growth conditions.
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
%   Synthesis_predictions : struct containing, at minimum:
%     .Mu_values                      : copy of Mu_vector
%     .Pi_pred_values                 : predicted synthesis-rate Π over Mu_vector
%     .(additional internal variables): host/circuit state quantities used in SynTwin
%     .Relative_sensitivities         : struct with relative sensitivities vs. parameters
%         .S_Pi_NA_values
%         .S_Pi_Omega_values
%         .S_Pi_RBS_k0_sigma0_values
%         .S_Pi_RBS_inv_sigma0_values
%
%   (Field names are kept consistent with downstream Results_Tensor construction
%    scripts shipped with SynTwin.)
%
% DEPENDENCIES
%   - SynTwin initialization is recommended before use:
%       ROOT = init_SynTwin(...);
%   - Requires the HEM surrogate struct (HEM) generated/provided with SynTwin.
%   - Used by estimation/result-generation workflows under Estimation_Pi/*.
%
% USAGE
%   S = Get_synthesis_predictions(HEM, model_c, Mu_vec, ...
%         RBS_k0_sigma0, RBS_inv_sigma0, Omega, Gen_cn);
%   Pi_pred = S.Pi_pred_values;
%   relS    = S.Relative_sensitivities;
%
% NOTES
%   - For optimization loops where only Π is needed, use the lightweight
%     Get_synthesis_predictions_lite.m for substantial speedups.
%   - Sensitivities returned here are relative (dimensionless) and are often
%     mapped to absolute sensitivities by multiplying by Π, as done in
%     Results-tensor generation scripts.
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
        Synthesis_predictions.Pi_pred_values(r,1) = Current_TU.Pi_pred_common_avg;
        Synthesis_predictions.fs_pred_values(r,1) = Current_TU.fs_common_avg;
        Synthesis_predictions.Mu_values(r,1) = Mu_vector(r);
        Synthesis_predictions.varphi_pred_values(r,1) = Current_TU.varphi_common_avg;
        Synthesis_predictions.load_pred_values(r,1) = Current_TU.theta_common_avg;
        % Estimation of the effective translation rate  KA_t_pred  = Pi_pred*dmA/N_A/omega_A:
        Synthesis_predictions.KA_t_pred_values(r,1)  = Synthesis_predictions.Pi_pred_values(r) *model_c.dm_c/Gen_cn/Omega;
        Synthesis_predictions.KRBS_c_values(r,1)  = KRBS_c(Synthesis_predictions.fs_pred_values(r),RBS_k0_sigma0,RBS_inv_sigma0);
       % Estimation of the relative (logarithmic) sensitivity of  the effective translation rate w.r.t. variations in the flux of resources:
       % Obtained as Svarphi_KA_t = 1/(1 + KRBS_c(f_s, k0, sigma0)/dmA*varphi)/varphi
       Synthesis_predictions.Svarphi_KA_t_values(r,1) = 1./(1+Synthesis_predictions.KRBS_c_values(r)./model_c.dm_c.*Synthesis_predictions.varphi_pred_values(r)  )./Synthesis_predictions.varphi_pred_values(r) ;
       % Estimation of the relative sensitivity of the synthesis rate w.r.t. biopart parameters:
       Synthesis_predictions.Relative_sensitivities.S_Pi_NA_values(r,1) = 1./Gen_cn;
       Synthesis_predictions.Relative_sensitivities.S_Pi_Omega_values(r,1) = 1./Omega;
       factor_sensitivity = 1./(RBS_inv_sigma0 + Current_TU.fs_common_avg + RBS_k0_sigma0*Current_TU.varphi_common_avg/model_c.dm_c);
       Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_k0_sigma0_values(r,1) = (RBS_inv_sigma0 + Current_TU.fs_common_avg)/RBS_k0_sigma0*factor_sensitivity;
       Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_inv_sigma0_values(r,1) = -factor_sensitivity;

end  % for values in Mu_exp_vector

end