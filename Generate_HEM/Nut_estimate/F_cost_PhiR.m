function J = F_cost_PhiR(x)
global Mu_exp_PhiR PhiR_exp
   
    predicted_PhiR =  x(1) + x(2).*Mu_exp_PhiR;

    sqr_prediction_error = ( predicted_PhiR -PhiR_exp ).^2  ;
    J = sum(sqr_prediction_error)/length(Mu_exp_PhiR);


 %  log_prediction_error = abs(log10(Nu_t_exp).*( log10(predicted_nut) -log10(Nu_t_exp) )) ;
  %  log_prediction_error = abs(abs(Nu_t_exp).^3.*( log10(predicted_nut) -log10(Nu_t_exp) )) ;
 %  abs_prediction_error = abs( predicted_nut -Nu_t_exp )  ;
  %  J = sum(log_prediction_error)/length(Mu_exp_nut);
   %  J = sum(abs_prediction_error)/length(Mu_exp_nut);
end