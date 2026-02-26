% eval_Host_fs_theta_v2
% SynTwin core routine: evaluate HEM steady-state predictions as a function
% of substrate profile (f_s) and HEM parameter vector (theta).
%
% DESCRIPTION
%   Given a substrate-dependent profile f_s and a parameter vector theta,
%   this function evaluates the Host Equivalent Model (HEM) steady-state
%   quantities required by the HEM fitting and surrogate construction
%   pipeline (e.g., predicted growth-related and allocation variables).
%
%   This function is used during HEM parameter fitting (WT/host calibration)
%   and during surrogate generation.
%
% INPUTS
%   model_p : struct
%       Host/model configuration parameters (constants and coefficients).
%   f_s     : double or vector
%       Substrate profile value(s) used to evaluate the host state.
%   theta   : double vector
%       HEM parameter vector being evaluated.
%
% OUTPUTS
%   Eval_Host_theta_results : struct
%       Struct containing the HEM predicted variables at the evaluated f_s.
%       (Fields depend on the current SynTwin HEM formulation.)
%
% USED BY
%   - Obtain_HEM_v1 (HEM surrogate fitting / generation)
%
% NOTES
%   - This routine is part of the HEM surrogate toolchain located under:
%       Generate_HEM/HEM_Surrogate/
%   - Typical SynTwin users do not need to call this function directly.
%
% SEE ALSO
%   Obtain_HEM_v1, Find_varphiHEM_thetaExo, Get_synthesis_predictions
function Eval_Host_theta_results = eval_Host_fs_theta_v2(model_p, f_s,theta) 

%Protein mass model 
mh = @(mu) ( model_p.c0 +  model_p.c1*mu.^model_p.c3)./( 1 + model_p.c2*mu.^model_p.c3 );

% nu_t model 
 nu_t_fs = @(f_s) ( model_p.nu_max*model_p.gamma1_fs*f_s./(model_p.gamma2_fs + f_s ) ); %aa/min

% RBS strengths
KRBS_r = @(f_s) model_p.k0_r./(1 + model_p.sigma0_r/model_p.nu_max*nu_t_fs(f_s) );     
KRBS_nr = @(f_s) model_p.k0_nr./(1 + model_p.sigma0_nr/model_p.nu_max*nu_t_fs(f_s) );   

 % Resources recruitment strengths
 J_r = @(f_s,varphy) model_p.Em_r*model_p.Omega_r/model_p.dm_r*KRBS_r(f_s)./(1 + KRBS_r(f_s)/model_p.dm_r.*varphy);  
 J_nr = @(f_s,varphy) model_p.Em_nr*model_p.Omega_nr/model_p.dm_nr*KRBS_nr(f_s)./(1 + KRBS_nr(f_s)/model_p.dm_nr.*varphy);     

mp_estimated =  mh(f_s*model_p.mu_max);
varphi_hat = model_p.varphi_profile; % Initial estimated profile to find the fixed point below
iter_Jr = NaN*ones(1,length(f_s)); %Preallocate space
iter_Jnr = NaN*ones(1,length(f_s)); %Preallocate space
iter_WSum_h = NaN*ones(1,length(f_s)); %Preallocate space
iter_Phi_r_w = NaN*ones(1,length(f_s)); %Preallocate space
iter_factor = NaN*ones(1,length(f_s)); %Preallocate space
varphi_hat_updated = NaN*ones(1,length(f_s)); %Preallocate space
mu_estimated = NaN*ones(1,length(f_s)); %Preallocate space
nut_vector= nu_t_fs(f_s);

for i=1:length(f_s)
  mp_fs_i =  mp_estimated(i);  %Estimated mass for the theoretical growth rate associated to the current substrate f_s value
  while true 
        % Estimate hat(varphi)_k = f( varphi_{k-1}, fs(i), model_params )
       iter_Jr(i) = J_r(f_s(i),varphi_hat(i) );
       iter_Jnr(i)  = J_nr(f_s(i),varphi_hat(i) );
       iter_WSum_h(i) = model_p.WEm_r *model_p.N_r *iter_Jr(i) + model_p.WEm_nr *model_p.N_nr*iter_Jnr(i) + theta;
       iter_Phi_r_w(i) = model_p.N_r *iter_Jr(i)./(1+iter_WSum_h(i)); 
       iter_factor(i)  = (model_p.Phi_m*iter_Phi_r_w(i))^2./(model_p.N_r *iter_Jr(i)  + model_p.N_nr*iter_Jnr(i) );
       varphi_hat(i)  = model_p.m_aa/model_p.ribosome_mass^2*iter_factor(i)*mp_fs_i.*nut_vector(i);
      % Use hat(varphi)_k to update hat(varphi)_updated = f( hat(varphi)_k, fs(i), model_params )
       iter_Jr(i) = J_r(f_s(i),varphi_hat(i) );
       iter_Jnr(i)  = J_nr(f_s(i),varphi_hat(i) );
       iter_WSum_h(i) = model_p.WEm_r *model_p.N_r *iter_Jr(i) + model_p.WEm_nr *model_p.N_nr*iter_Jnr(i) + theta;
       iter_Phi_r_w(i) = model_p.N_r *iter_Jr(i)./(1+iter_WSum_h(i)); 
       iter_factor(i)  = (model_p.Phi_m*iter_Phi_r_w(i))^2./(model_p.N_r *iter_Jr(i)  + model_p.N_nr*iter_Jnr(i) );
       varphi_hat_updated(i)  = model_p.m_aa/model_p.ribosome_mass^2*iter_factor(i)*mp_fs_i.*nut_vector(i);
       % Check for the difference hat(varphi)_updated - hat(varphi)_k
       if abs(varphi_hat_updated(i)  - varphi_hat(i)) < 1e-4
            % We have found the fixed point
            varphi_hat(i)  = varphi_hat_updated(i) ;
            mu_estimated(i)  = 1/(model_p.N_r*model_p.lp_r)*model_p.Phi_m*iter_Phi_r_w(i).*nut_vector(i);
            break; %Leave the while loop
       else
           % Iterate again
             varphi_hat(i)  = varphi_hat_updated(i) ; %Iterate with the new estimated varphi profile
       end
   end %while
end % for f_s(i)
% Return values 
Eval_Host_theta_results.varphi_estimated= varphi_hat;
Eval_Host_theta_results.mu_estimated = mu_estimated;
Eval_Host_theta_results.WSum_h = iter_WSum_h;
Eval_Host_theta_results.Sum_h = model_p.N_r *iter_Jr+ model_p.N_nr*iter_Jnr;
Eval_Host_theta_results.J_r = iter_Jr;
Eval_Host_theta_results.J_nr = iter_Jnr;
% free ribosomes r=varphi/mu
Eval_Host_theta_results.free_r = Eval_Host_theta_results.varphi_estimated./Eval_Host_theta_results.mu_estimated;
% available mature ribosomes   r_a = r*(1+ WSum_h) = varphi/mu*(1+ WSum_h)
Eval_Host_theta_results.available_r_a = (1+ Eval_Host_theta_results.Sum_h).*Eval_Host_theta_results.free_r ;
% active translating ribosomes r_t = Sum_h/(1+ WSum_h)*r_a = Sum_h*r
Eval_Host_theta_results.active_r_t = Eval_Host_theta_results.Sum_h.*Eval_Host_theta_results.free_r;
%  total ribosomes r_T = r_a/Phi_m 
Eval_Host_theta_results.total_r_T = Eval_Host_theta_results.available_r_a/model_p.Phi_m;
end

  



   
     
   
 
 
 