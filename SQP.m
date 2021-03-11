function [x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval, rho_all, nb_eval_all] = SQP(x, lambda, ...
                                problem, merite, choix, bornes, max_iter, max_eval, rho, eps, rho_max)
%%%---- Algorithme SQP ----%%%                            

%%%------------- Initialisation ---------------------------------------%%%

    % taille du vecteur d'entrée
    n = length(x); 
    % nombre d'iterations de SQP
    nb_iter = 0; 
    % nombre d'appels de 'problem'
    nb_eval = 0; 
    % pas pour le gradient par différences finies unilatérales
    h = repmat(1e-8, n, 1); 
    h = h.*x; 
    % initialisation de la hessienne
    H = eye(n); 
    % bornes des variables
    borne_inf = bornes(:,1);
    borne_sup = bornes(:,2);
    % indicateur pour le résultat de la globalisation
    outcome = 0; 
    % initialisation du stokage des valeurs par itération:
    x_all = []; 
    f_all = [];
    c_all = [];
    lambda_all = [];
    grad_L_norm_all = [];
    rho_all = [];
    nb_eval_all = [];
    
%%%------------------ Coeur de l'algo ---------------------------------%%%

    nb_iter = nb_iter + 1;
    %fprintf("Iteration n°: %d\n", nb_iter);    
    
    [f, c] = problem(x); % f = f(x); c = c(x)
    nb_eval = nb_eval + 1;
    % gradient par différences finies
    [gradf, gradc] = grad_diff_finies(x, h, problem);
    nb_eval = nb_eval + n + 1;
    % gradient du lagrangien
    grad_L = gradient_lagrangien(lambda, gradf, gradc);
    % critère d'arrêt        
    stop = (nb_iter >= max_iter) | (nb_eval >= max_eval) | (norm(grad_L) < eps);

    while ~stop
        nb_iter = nb_iter + 1;
        %fprintf("Iteration n°: %d\n", nb_iter);
        
        x_all = [x_all, x]; 
        lambda_all = [lambda_all, lambda];
        rho_all = [rho_all, rho];

        if (outcome == 0)
            if (nb_iter > 2)
               % à partir de la deuxième itération on approxime la
               % hessienne par quasi-newton, en utilisant BFGS ou SR1 
               H = quasi_newton(H, x, x_avant, grad_L, grad_L_avant, choix);
            end
        
            % modification hessien pour la rendre définie positive
            H = modif_hess(H, 0.001);
        
            % solution problème quadriatique
            [d_QP, lambda] = solution_quad(gradf, gradc, c, H);
     
            % on stocke les valeurs d'avant
            x_avant = x;
            f_avant = f;
            grad_L_avant = gradient_lagrangien(lambda, gradf, gradc);
        end
        
        % on stocke les valeurs de la fonction, la contrainte et la
        % norme du gradient du lagrangien 
        f_all = [f_all, f];
        c_all = [c_all, c];
        grad_L_norm_all = [grad_L_norm_all, norm(grad_L)];
        
        % globalisation, pour améliorer la descente
        [x, f, c, nb_eval_glob, outcome] = globalisation(x, problem, ...
            merite, f, c, gradf, d_QP, borne_inf, borne_sup, rho, outcome);
        nb_eval = nb_eval + nb_eval_glob;
        nb_eval_all = [nb_eval_all, nb_eval];
        
        % si échec de type 1 : dérivée directionnelle positive
        if outcome == 1
            H = eye(n);
            outcome = 2;
            %disp("hessienne réinitialisée");
            stop = (nb_iter >= max_iter) | (nb_eval >= max_eval);
            continue;
        end
        
        % si échec de type 2 : dérivée positive malgré réinitialisation de
        %                      la hessienne 
        if outcome == 2 && rho/2 <= rho_max
            rho = 2*rho; % si 2ème échec d'affilée on augmente rho
            %fprintf("On augmente rho!\nrho = %d\n", rho);
            %outcome = 1;
            stop = (nb_iter >= max_iter) | (nb_eval >= max_eval);
            continue;
        end
        
        if outcome == 0
            [gradf, gradc] = grad_diff_finies(x, h, problem);
            nb_eval = nb_eval + n + 1;
            grad_L = gradient_lagrangien(lambda, gradf, gradc);
        end
    
        % critères d'arrêt
        stop = (nb_iter >= max_iter) | (nb_eval >= max_eval) ...
                 | (norm(grad_L) < eps) | (norm(x-x_avant) < eps) ...
                 | (abs(f - f_avant) < eps);
    end
    
    % on garde les toutes dernières valeurs
    x_all = [x_all, x];
    lambda_all = [lambda_all, lambda];
    rho_all = [rho_all, rho];
    f_all = [f_all, f];
    c_all = [c_all, c];
    grad_L_norm_all = [grad_L_norm_all, norm(grad_L)];
    nb_eval_all = [nb_eval_all, nb_eval];

    %fprintf("norm(grad_L) = %f\n", norm(grad_L));
    %fprintf("norm(x - x_avant) = %f\n", norm(x-x_avant));
    %fprintf("norm(f - f_avant) = %f\n", norm(f - f_avant));
end