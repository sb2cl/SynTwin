function J = F_cost_nut(x,Mu_exp_nut,Nu_t_exp)
global Mu_exp_nut Nu_t_exp

   
    predicted_nut =  1e3*x(1).*Mu_exp_nut./( x(2) + Mu_exp_nut );
    sqr_prediction_error = ( predicted_nut -Nu_t_exp ).^2  ;
    J = sum(sqr_prediction_error)/length(Mu_exp_nut);
end