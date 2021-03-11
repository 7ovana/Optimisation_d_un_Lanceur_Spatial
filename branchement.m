function [f, c] = branchement(THETA, simulation_handle, R_c)
    
    [R_tf, V_tf] = simulation_handle(THETA);
    % fonction à minimiser
    f = - norm(V_tf);
    % contrainte 
    c = [norm(R_tf) - R_c; R_tf'*V_tf];
end