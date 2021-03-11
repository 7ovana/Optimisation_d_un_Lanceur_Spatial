function [f, c, bornes, Mi] = ariane(m, Vp, m_satelite, k, v_e)
%%%--- Notre cas : Ariane 6 ---%%%

    % masse sèche
    ms = k.*m;        
    
    Mi4 = m_satelite;  
    
    Mf3 = Mi4 + ms(3);
    Mi3 = Mf3 + m(3);
    
    Mf2 = Mi3 + ms(2);
    Mi2 = Mf2 + m(2);
    
    Mf1 = Mi2 + ms(1);
    Mi1 = Mf1 + m(1);
    
    Mij = [Mi1; Mi2; Mi3];
    Mfj = [Mf1; Mf2; Mf3];
    
    V = sum(v_e.*log(Mij ./ Mfj));
    
% Fonction à minimiser, f
    f = Mi1;
    
% Contrainte d'égalité, c    
    c = V - Vp;
    
    if nargout > 2
        bornes = [8000, 50000; 5000, 20000; 4000, 10000];
    end
    
    if nargout > 3
        Mi = Mij;
    end
end
