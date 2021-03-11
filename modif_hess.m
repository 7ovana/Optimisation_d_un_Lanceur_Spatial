function H_modif = modif_hess(H,eps)
%Rendre H SDP si elle ne l'est pas
%valeurs propres
Vp = eig(H); 
I = eye(length(Vp));
minVp = min(Vp);
% VP strictement positives on ne fait rien
if (minVp > eps)
    H_modif = H;
else 
    %VP negatives on ajoute la valeur absolue du minimum des VP + eps pour
    %la rendre SDP
    tau = abs(minVp) + eps;
    H_modif = H + tau*I ;
end
end 