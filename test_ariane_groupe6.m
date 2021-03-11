function [x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval, rho_all] = test_ariane_groupe6(depart_newton, V_p, choix, max_iter, max_eval, rho, eps, max_rho)
%%%---- TEST: Comparaison des valeurs obtenues analytiquement avec les 
%%%           valeurs obtenues par SQP, pour les masses des ergols 
%%%           pour nos données du problème Ariane

%%% ---- DONNEES DU PROBLEME ---- %%%
    
    % indice constructif par étage
    k = [0.10; 0.15; 0.20];
    % vitesse d'ejection par étage (m/s)
    v_e = [2600; 3000; 4400];
    % masse du satelite
    m_satelite = 2000;
    
% Solution analytique
    sol_anal = ariane_solution_analytique(depart_newton, k, v_e, m_satelite, V_p);
    [f_anal, c_anal, bornes] = ariane_groupe6(sol_anal, V_p, m_satelite, k, v_e);
    
% Initialisation pour l'optimisateur SQP 
    x_init = [15000; 13000; 5000];
    lambda_init = 0;

    % Rendre ariane_groupe6 fonction d'une seule variable x
    ariane6 = @(x) ariane_groupe6(x, V_p, m_satelite, k, v_e);

% Appel à SQP
    [x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval, rho_all] = SQP(x_init, ...
        lambda_init, ariane6, @merite, choix, bornes, max_iter, max_eval, rho, eps, max_rho);
 
% Résultat de SQP
    sol = x_all(:,end);
    [f_sol, c_sol] = ariane6(sol);
       
    fprintf("nb_iter = %d\nnb_eval = %d\n\n", nb_iter, nb_eval);
    fprintf("rho max atteint?  rho = %d\n\n", rho_all(end));

    fprintf("sol_SQP = [%f; %f; %f]\n", sol);
    fprintf("f(sol_SQP) = %f\n", f_sol);
    fprintf("c(sol_SQP) = %f\n\n", c_sol);
    
    fprintf("sol_analytique = [%f; %f; %f]\n", sol_anal);
    fprintf("f(sol_analytique) = %f\n", f_anal);
    fprintf("c(sol_analytique) = %f\n\n", c_anal);
    
    fprintf("f_SQP - f_analytique = %f\n", f_sol - f_anal);
    fprintf("c_SQP - c_analytique = %f\n", abs(c_sol - c_anal));
    fprintf("sol_SQP - sol_analytique = [%f; %f; %f]\nn", abs(sol - sol_anal));
end