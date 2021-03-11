function Hk = quasi_newton(H_0, x, x0, grad_L, grad_L0, choix)
%Renvoie une approximation de Hk à partir du dernier deplacement
%H_0 correspond H_k-1
%x correspond à x_k et x0 à x_k-1
%grad_l gradient du lagrangien au point x
%grad_l0 gradient du lagrangien au point x0
% choix = 'SR1' pour SR1 ou 'BFGS' pour BFGS

%variation de x
d = x - x0;
y = grad_L - grad_L0;

% Première méthode : SR1
if (choix == "SR1")
    if (d'*(y - H_0*d) ~= 0)
        Hk = H_0 + ((y-H_0*d)*(y-H_0*d)')/(d'*(y-H_0*d));
    else
        Hk = H_0;
    end
% Deuxième méthode : BFGS    
elseif (choix == "BFGS")
   if (y'*d > 0)
        Hk = H_0 + (y*y')/(y'*d) - (((H_0*d)*d') *H_0)/(d'*H_0*d);
   else
        Hk = H_0;
   end
else
    return;
end
end


