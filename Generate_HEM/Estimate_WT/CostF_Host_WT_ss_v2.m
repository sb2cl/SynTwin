function F = CostF_Host_WT_ss_v2(x)
%% CostF_Host_WT_ss_v2
% Objective (cost) function for steady-state WT-host parameter estimation.
%
% Used by: Estim_Runs_HostWT_v2 (MEIGO/eSS)
%
% Given a candidate parameter vector x, this function evaluates the mismatch
% between WT-host model predictions and the collated literature datasets
% packaged with SynTwin (see input_syn_data.m).
%
% INPUT
%   x : parameter vector to be estimated by the optimizer
%
% OUTPUT
%   F : scalar objective value (lower is better)


global model_p;

k0_r = x(1,1);
k0_nr = x(1,2);
sigma0_r = x(1,3);
sigma0_nr = x(1,4);
Omega_r = x(1,5);
Omega_nr = x(1,6);
Phi_m= x(1,7);

%Protein mass model m = ( 115.355e-15 +  455.799e-9*μ.^3.3448)./( 1 + 1028.37e3.*μ^3.3448 ) 
mh = @(mu) ( model_p.c0 +  model_p.c1*mu.^model_p.c3)./( 1 + model_p.c2*mu.^model_p.c3 );

% nu_t model 
 nu_t_fs = @(f_s) ( model_p.nu_max*model_p.gamma1_fs*f_s./(model_p.gamma2_fs + f_s ) ); %aa/min

% RBS strengths
KRBS_r = @(f_s) k0_r./(1 + sigma0_r/model_p.nu_max*nu_t_fs(f_s) );     
KRBS_nr = @(f_s) k0_nr./(1 + sigma0_nr/model_p.nu_max*nu_t_fs(f_s) );   

 % Resources recruitment strengths
 J_r = @(f_s,varphy) model_p.Em_r*Omega_r/model_p.dm_r*KRBS_r(f_s)./(1 + KRBS_r(f_s)/model_p.dm_r.*varphy);  
 J_nr = @(f_s,varphy) model_p.Em_nr*Omega_nr/model_p.dm_nr*KRBS_nr(f_s)./(1 + KRBS_nr(f_s)/model_p.dm_nr.*varphy);     

f_s = model_p.f_s;
mu_exp = model_p.mu_test;

mp_estimated =  mh(f_s*model_p.mumax);
varphi_hat = model_p.varphi_profile; % Initial estimated profile to find the fixed point below
iter_Jr = NaN*ones(1,length(f_s)); %Preallocate space
iter_Jnr = NaN*ones(1,length(f_s)); %Preallocate space
iter_WSum_h = NaN*ones(1,length(f_s)); %Preallocate space
iter_Phi_t = NaN*ones(1,length(f_s)); %Preallocate space
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
       iter_WSum_h(i) = model_p.WEm_r *model_p.N_r *iter_Jr(i) + model_p.WEm_nr *model_p.N_nr*iter_Jnr(i);
       iter_Phi_t(i) = model_p.N_r *iter_Jr(i)./(1+iter_WSum_h(i)); 
       iter_factor(i)  = (Phi_m*iter_Phi_t(i))^2./(model_p.N_r *iter_Jr(i)  + model_p.N_nr*iter_Jnr(i) );
       varphi_hat(i)  = model_p.m_aa/model_p.ribosome_mass^2*iter_factor(i)*mp_fs_i.*nut_vector(i);
      % Use hat(varphi)_k to update hat(varphi)_updated = f( hat(varphi)_k, fs(i), model_params )
       iter_Jr(i) = J_r(f_s(i),varphi_hat(i) );
       iter_Jnr(i)  = J_nr(f_s(i),varphi_hat(i) );
       iter_WSum_h(i) = model_p.WEm_r *model_p.N_r *iter_Jr(i) + model_p.WEm_nr *model_p.N_nr*iter_Jnr(i);
       iter_Phi_t(i) = model_p.N_r *iter_Jr(i)./(1+iter_WSum_h(i)); 
       iter_factor(i)  = (Phi_m*iter_Phi_t(i))^2./(model_p.N_r *iter_Jr(i)  + model_p.N_nr*iter_Jnr(i) );
       varphi_hat_updated(i)  = model_p.m_aa/model_p.ribosome_mass^2*iter_factor(i)*mp_fs_i.*nut_vector(i);
       % Check for the difference hat(varphi)_updated - hat(varphi)_k
       if abs(varphi_hat_updated(i)  - varphi_hat(i)) < 1e-4
            % We have found the fixed point
            varphi_hat(i)  = varphi_hat_updated(i) ;
            mu_estimated(i)  = 1/(model_p.N_r*model_p.lp_r)*Phi_m*iter_Phi_t(i).*nut_vector(i);
            break; %Leave the while loop
       else
           % Iterate again
             varphi_hat(i)  = varphi_hat_updated(i) ; %Iterate with the new estimated varphi profile
       end
   end %while
end % for f_s(i)

% The cost function to be optimized is the sum of the absolute prediction
% errors between the experimental and the estimated growth rates:

F = sum(abs(mu_estimated - mu_exp).^2);


end

  



   
     
   
 
 
 
