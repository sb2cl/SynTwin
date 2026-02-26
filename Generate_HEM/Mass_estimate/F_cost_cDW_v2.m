function J = F_cost_cDW_v2(parameters,Mu_vector,Exp_cDW_vector)
    % MODEL  (Hill)
    predicted_mh = ( parameters(1) + 1e6*parameters(2).*Mu_vector.^parameters(4) )./( 1 + 1000*parameters(3).*Mu_vector.^parameters(4) );
    sq_prediction_error = (predicted_mh-Exp_cDW_vector).^2;
    J = sum(sq_prediction_error)/length(Mu_vector);
end