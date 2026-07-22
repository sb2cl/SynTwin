function J = F_cost_PhiR(x)
global Mu_exp_PhiR PhiR_exp
   
    predicted_PhiR =  x(1) + x(2).*Mu_exp_PhiR;

    sqr_prediction_error = ( predicted_PhiR -PhiR_exp ).^2  ;
    J = sum(sqr_prediction_error)/length(Mu_exp_PhiR);

end