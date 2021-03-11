function F = merite(f, c, rho)
%%%--- Fonction mérite ---%%%

%f : fonction a minimiser
%c : fonction contrainte
%rho : penalisation

F = f + rho*norm(c,1);
end