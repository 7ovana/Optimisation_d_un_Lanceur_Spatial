function [x_opt, f_next, c_next, nb_eval_glob, outcome] = globalisation(x, problem, ...
                             merite, f, c, gradf, d_QP, borne_inf, borne_sup, rho, outcome)
 %%%---- GLOBALISATION ----%%%
 
 %%%------------- Initialisation -------------------------------------------------%%%

    % pas d'Armijo
    s = 1; 
    % condition de décroissance suffisante
    c1 = 0.1; 
    % nombre d'appels à la fonction 'problem'
    nb_eval_glob = 0; 
    % fonction de mérite en le x courant       
    F_x = merite(f, c, rho); 
    nb_eval_glob = nb_eval_glob + 1; 
    % indicateur sur la réussite de la recherche linéaire d'Armijo
    armijo = false; 
    % initialisation du x optimal comme le x courant
    x_opt = x;  
    % initialisation de f_next et c_next comme le f et c courant
    f_next = f; c_next = c;
    % calcul de la derivée directionnelle
    derivee_direction = gradf'*d_QP - rho*norm(c,1); 
    
%%%------------------ Coeur de l'algo -------------------------------------------%%%

    if derivee_direction >= 0 % si la derivee est positive
        
        if outcome == 2 % si on avait déjà réinitialisé la hessienne
            %fprintf("outcome = %d, hessienne déjà réinitialisée\n", outcome);
            return; % on sort pour augmenter rho
        end
        
        outcome = 1; % sinon on note l'échec en mettant outcome à 1
        %fprintf("outcome = %d, derivee positive\n", outcome); 
        return; % on sort pour réinitialiser la hessienne
    
    else % si la derivée est bien négative, on fait la recherche linéaire
        while (s > 1e-6)         
            xsd = projection_bornes(x + s*d_QP, borne_inf, borne_sup);
            [fsd, csd] = problem(xsd);
            nb_eval_glob = nb_eval_glob + 1;
            F_xsd = merite(fsd, csd, rho);
            
            if (F_xsd < F_x + c1*s*derivee_direction) 
                armijo = true;
                break;
            end
            % on divise s par deux tant qu'on trouve le bon
            s = s/2;          
        end
         
         if armijo % si on réussit à trouver un bon pas d'Armijo
            outcome = 0;
            x_opt = xsd; % on garde ce xsd pour le pas prochain de SQP
            f_next = fsd; c_next = csd; % et on garde la valeur de f et c en xsd 
            %fprintf("Armijo success: outcome = %d,  s = %f\n", outcome, s);
            return;
         else 
            if outcome == 0 % si la dernière fois qu'on était ici, on avait réussit
                outcome = 1; % maintenant on marque l'échec car s est trop petit
            end              % si outcome est déjà 2, on le laisse ainsi pour augmenter rho
            %fprintf("s trop petit, on sort de globalisation, outcome = %d\n", outcome);
            return;
         end
     end
end