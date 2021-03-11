function [y, f_prim] = equation_analytique(x, k, v, Vp)
%%%--- Equation pour la résolution analytique du problème d'étagement ---%%%   
    Omega = k ./ (1 + k);
    d1 = v(1) - v(3);
    d2 = v(2) - v(3);
    eta = v .* Omega; 
    
% l'équation, vu comme une fonction de x:
    y = v(1)*log(eta(3)*x + d1) + v(2)*log(eta(3)*x + d2) + v(3)*log(x) ...
       - v(1)*log(eta(1)) - v(2)*log(eta(2)) - Vp;
   
% la dérivée de cette fonction de x:
    f_prim = v(1)*eta(3) / (eta(3)*x + d1) + v(2)*eta(3) / (eta(3)*x + d2) + v(3) / x;
 
end
