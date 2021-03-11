function [x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval, rho_all, nb_eval_all] = test_MHW4D(choix, max_iter, max_eval, rho, eps, max_rho)
%%% --- Test du cas test 1 : MHW4D ---%%%

% Solution donnée dans le polycopié
    sol_poly = [-1.2366; 2.4616; 1.1911; -0.2144; -1.6165];
    [f_poly, c_poly, bornes] = MHW4D(sol_poly);
    disp('bornes pour x:');
    disp(bornes);
    
% Initialisation
    x_init = [-1; 2; 1; -2; -2];
    lambda_init = [0; 0; 0];
    
% Appel à SQP
    [x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval, rho_all, nb_eval_all] = SQP(x_init, lambda_init, ...
                                @MHW4D, @merite, choix, bornes, max_iter, max_eval, rho, eps, max_rho);

% Notre solution                            
    sol = x_all(:,end);
    f_sol = f_all(:,end);
    c_sol = c_all(:,end);
    
    fprintf("nb_iter = %d\nnb_eval = %d\n\n", nb_iter, nb_eval);
    fprintf("rho max atteint?  rho = %d\n\n", rho_all(end));

    fprintf("sol = [%f; %f; %f; %f; %f]\n", sol);
    fprintf("f(sol) = %f\n", f_sol);
    fprintf("c(sol) = [%f; %f; %f]\n", c_sol);
    
    fprintf("sol_poly = [%f; %f; %f; %f; %f]\n", sol_poly);
    fprintf("f(sol_poly) = %f\n", f_poly);
    fprintf("c(sol_poly) = [%f; %f; %f]\n\n", c_poly);
    
    fprintf("sol - sol_poly =  [%f; %f; %f; %f; %f]\n", sol - sol_poly);
    fprintf("f - f_poly = %f\n", f_sol - f_poly);
    fprintf("c - c_poly = [%f; %f; %f]\n", c_sol - c_poly);
end