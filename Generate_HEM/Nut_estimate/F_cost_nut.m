function J = F_cost_nut(x,Mu_exp_nut,Nu_t_exp)
global Mu_exp_nut Nu_t_exp

   
    predicted_nut =  1e3*x(1).*Mu_exp_nut./( x(2) + Mu_exp_nut );
 %  log_prediction_error = abs(log10(Nu_t_exp).*( log10(predicted_nut) -log10(Nu_t_exp) )) ;
  %  log_prediction_error = abs(abs(Nu_t_exp).^3.*( log10(predicted_nut) -log10(Nu_t_exp) )) ;
    sqr_prediction_error = ( predicted_nut -Nu_t_exp ).^2  ;
 %  abs_prediction_error = abs( predicted_nut -Nu_t_exp )  ;
  %  J = sum(log_prediction_error)/length(Mu_exp_nut);
    J = sum(sqr_prediction_error)/length(Mu_exp_nut);
   %  J = sum(abs_prediction_error)/length(Mu_exp_nut);
end