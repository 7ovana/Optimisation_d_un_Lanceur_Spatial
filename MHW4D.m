function [f, c, bornes] = MHW4D(x)
%%%--- Cas test 1 : MHW4D ---%%%

% Fonction à minimiser, f
    f = (x(1) - 1)^2 + (x(1) - x(2))^2 + (x(2) - x(3))^3 + ...
        (x(3) - x(4))^4 + (x(4) - x(5))^4;
    
% Contrainte d'égalité, c
    c(1,:) = x(1) + x(2)^2 + x(3)^2 - 3*sqrt(2) - 2;
    c(2,:) = x(2) - x(3)^2 + x(4) - 2*sqrt(2) + 2;
    c(3,:) = x(1)*x(5) - 2;    
    
% si 3 arguments en sortie, on donne les bornes des variables
    if nargout == 3
        bornes = [-2, 0; 1, 3; 0, 2; -3, 0; -3, -1];
    end
end
    