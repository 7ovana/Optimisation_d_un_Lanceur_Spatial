%%% REMARQUES

% Appels utiles:

% MH4WD

[x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval, rho_all] = test_MHW4D('SR1', 50, 1e4, 6, 1e-5, 2^20);
[x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval, rho_all] = test_MHW4D('BFGS', 50, 1e4, 6, 1e-5, 2^20);

% Arianne 1

[x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval] = test_ariane1('SR1', 1e3, 1e7, 8, 1e-4, 2^20);
[x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval] = test_ariane1('BFGS', 1e3, 1e7, 8, 1e-4, 2^20);

% Arianne 6

% V_p
[x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval] = test_ariane_groupe6(4, V_p, 'SR1', 1e3, 1e7, 2, 1e-4, 2^12);
[x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval] = test_ariane_groupe6(4, V_p, 'BFGS', 1e3, 1e7, 2, 1e-4, 2^12);

% V_p + deltaV
[x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval] = test_ariane_groupe6(4, V_p + deltaV, 'SR1', 1e3, 1e7, 2, 1e-4, 2^12);
[x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval] = test_ariane_groupe6(4, V_p + deltaV, 'BFGS', 1e3, 1e7, 2, 1e-4, 2^12);





