function [f, c, bornes, Mi] = ariane1(m)
%%%--- Cas test 2 : Ariane 1 ---%%%

    % indice constructif par étage
    k = [0.1101; 0.1532; 0.2154]; 
    % vitesse d'éjection par étage
    ve = [2647.2; 2922.4; 4344.3];
    % masse sèche
    ms = k.*m;        
    % masse du satelite
    mu = 1700;
    % vitesse propulsive réquise
    deltaV_requis = 11527;

    Mi4 = mu;  
    
    Mf3 = Mi4 + ms(3);
    Mi3 = Mf3 + m(3);
    
    Mf2 = Mi3 + ms(2);
    Mi2 = Mf2 + m(2);
    
    Mf1 = Mi2 + ms(1);
    Mi1 = Mf1 + m(1);
    
    Mij = [Mi1; Mi2; Mi3];
    Mfj = [Mf1; Mf2; Mf3];
    
    V = sum(ve.*log(Mij ./ Mfj));
    
% Fonction à minimiser, f
    f = Mi1;
    
% Contrainte d'égalité, c    
    c = V - deltaV_requis;
    
    if nargout > 2
        bornes = [100000, 150000; 20000, 50000; 5000, 10000];
    end
    
    if nargout > 3
        Mi = Mij;
    end
end
