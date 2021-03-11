function nablah_L = gradient_lagrangien(lambda, gradf, gradc)
% Gradient du lagrangien
    nablah_L = gradf + gradc*lambda;
end