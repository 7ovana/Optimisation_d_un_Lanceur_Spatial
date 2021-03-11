function [x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval, rho_all, nb_eval_all] = test_ariane1(choix, max_iter, max_eval, rho, eps, max_rho)
%%%--- Test de notre algo SQP ---%%%

% Solution donnée dans le polycopié
    sol_ariane = [145349; 31215; 7933];
    [f_ariane, c_ariane, bornes] = ariane1(sol_ariane);
    disp('bornes pour x:');
    disp(bornes);
    
% Initialisation 
    %x_init = [140000; 25000; 6500];
    x_init = [130000; 30000; 9000];
    lambda_init = 0;

% Appel à SQP
    [x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval, rho_all, nb_eval_all] = SQP(x_init, ...
        lambda_init, @ariane1, @merite, choix, bornes, max_iter, max_eval, rho, eps, max_rho);

% Notre solution
    sol = x_all(:,end);
    f_sol = f_all(:,end);
    c_sol = c_all(:,end);
    
    fprintf("nb_iter = %d\nnb_eval = %d\n\n", nb_iter, nb_eval);
    fprintf("rho max atteint?  rho = %d\n\n", rho_all(end));

    fprintf("sol = [%f; %f; %f]\n", sol);
    fprintf("f(sol) = %f\n", f_sol);
    fprintf("c(sol) = %f\n\n", c_sol);
    
    fprintf("sol_ariane = [%f; %f; %f]\n", sol_ariane);
    fprintf("f(sol_ariane) = %f\n", f_ariane);
    fprintf("c(sol_ariane) = %f\n\n", c_ariane);
    
    fprintf("f - f_ariane = %f\n", f_sol - f_ariane);
    fprintf("c - c_ariane = %f\n", abs(c_sol - c_ariane));
    fprintf("sol - sol_ariane = [%f; %f; %f]\nn", abs(sol - sol_ariane));
    
end
