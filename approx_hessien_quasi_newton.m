function Q = approx_hessien_quasi_newton(x, h, problem)
    [f, c] = problem(x);
    lambda = 0;
    [gradf, gradc] = grad_jac_diff_finies(x, h, problem);
    grad_lagrangien = [gradf + gradc*lambda; c];
    
    % initialisation:
    H = eye(length(x));
    
    it = 0;
    while() % reinitialisation toutes les n iterations
        %BFGS:
        x_avant = x;
        lambda_avant = lambda;
        % x = faut updater x;
        % lambda = faut updater lambda;
        if (y'*d > 0)
            H = H + (y*y') / (y'*d) - (H*d*d'*H) / (d'*H*d);
        end
        % x = x + d_x;
        
        % est-ce qu'on utilise solution_quad pour lambda et d ???
        
        lambda = lambda + d_lambda;
        y = grad_lagrangien(x,lambda) - grad_lagrangien(x_avant, lambda_avant);
        % d = x - x_avant;
        it = it+1;
        
        if (it == n) 
            H = eye(length(x));
        end
    end
    