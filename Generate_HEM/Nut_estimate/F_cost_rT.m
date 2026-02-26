function J = F_cost_rT(x,Mu_vector,rT_new)
global rT_new
global Mu_vector
global model_num    
 if model_num==1
    predicted_rT =  1e5*x(1).*Mu_vector.^x(2);
elseif model_num==2
    predicted_rT = ( 1e3*x(1) + 1e9*x(2).*Mu_vector.^x(4) )./( 1 + 1e4*x(3).*Mu_vector.^x(4) );
end
    sq_prediction_error = (predicted_rT-rT_new).^2;
    J = sum(sq_prediction_error)/length(Mu_vector);
end