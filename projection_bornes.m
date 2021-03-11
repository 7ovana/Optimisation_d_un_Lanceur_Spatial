function x = projection_bornes(x, borne_inf, borne_sup)
%--- Projection dans les bornes autorisées pour x ---%
n = length(x);
for i = 1:n
    if(x(i) <= borne_inf(i)) % si en dessous du minimum autorisé
        x(i) = borne_inf(i); % on le projette sur le minimum
    elseif (x(i) >= borne_sup(i)) % si au-dessus du maximum autorisé
        x(i) =  borne_sup(i); % on le projette sur le maximum
    end
end
end

    