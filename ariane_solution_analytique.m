function [sol_anal, X] = ariane_solution_analytique(depart_newton, k, v_e, m_satelite, V_p)
%%% ----- SOLUTION ANALYTIQUE DU PROBLEME D'ETAGEMENT ---- %%%
%
% depart_newton: point de départ pour la méthode de Newton
% k: indice constructif par étage
% v_e: vitesse d'éjection par étage
% m_satelite: masse du satelite
% V_p: vitesse propulsive

    Omega = k ./ (1 + k);
    d1 = v_e(1) - v_e(3);
    d2 = v_e(2) - v_e(3);
    eta = v_e .* Omega; 
    Mi4 = m_satelite;
    
    % Ecriture des x_j
    equation = @(x) equation_analytique(x, k, v_e, V_p);
    [x3, ~] = newton(depart_newton, equation, 1e-8);
    
    fprintf("x_optimal = %f\n\n", x3);
    x2 =  1 / eta(2) * (eta(3)*x3 + d2);
    x1 =  1 / eta(1) * (eta(3)*x3 + d1);

    %Ecriture des m_ej
    Mi3 = (Mi4 * x3) / (1 + k(3) - k(3)*x3);
    Mi2 = (Mi3 * x2) / (1 + k(2) - k(2)*x2);
    Mi1 = (Mi2 * x1) / (1 + k(1) - k(1)*x1);

    Mf1 = Mi1 / x1;
    Mf2 = Mi2 / x2;
    Mf3 = Mi3 / x3;

    m_e1 = Mi1 - Mf1;
    m_e2 = Mi2 - Mf2;
    m_e3 = Mi3 - Mf3;
    
    sol_anal = [m_e1; m_e2; m_e3];
    X = [x1; x2; x3];
end