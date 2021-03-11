%%% --------------------- LOGICIEL FINAL -------------------------------%%%

%%%------------------- CONSTANTES UTILISEES ----------------------------%%%

% constante gravitationnelle terrestre (m^3/s^2)
mu = 3.986*1e14;
% coefficient de trainée du lanceur
c_x = 0.1;
% rayon terrestre 
R_terre = 6378137;

%%%----------------------- DONNEES DU PROBLEME -------------------------%%%

% hauteur cible
H_c = 250000;
% rayon de l'orbite
R_c = R_terre + H_c;
% vitesse cible sur l'orbite
V_c = sqrt(mu / R_c);
% masse du satelite 
m_satellite = 2000;
% indice constructif par étage
k = [0.10; 0.15; 0.20];
% vitesse d'ejection par étage (m/s)
v_e = [2600; 3000; 4400];
%acceleration initiale
alpha = [15; 10; 10]; 

%%%-------------------------- INITIALISATION ---------------------------%%%

% vitesse propulsive
V_p = V_c;
% écart entre la vitesse cible et la vitesse réelle
deltaV = 0.2*V_c;
% nombre d'itérations
nb_iter_total = 0;

% ------ ETAGEMENT ---------- % 

% bon x de départ pour SQP 
m_e_init = [15000; 13000; 5000];
% masse sèche initiale
m_s = k.*m_e_init;
% initialisation de lambda
lambda_m_e_init = 0;
% bornes des masses des ergols
bornes_m_e = [8000, 50000; 5000, 20000; 4000, 10000];
[M, ~, ~, Mi] = ariane(m_e_init, V_p + deltaV, m_satellite, k, v_e);

% ------ TRAJECTOIRE -------- %

% bon theta de départ pour SQP
theta_0 = 4.5; theta = [2; 6.5; 4];
THETA_init = [theta_0; theta];
% initialisation de lambda
lambda_THETA_init = [0; 0];
% bornes de THETA
bornes_THETA = [-90, 90; -90, 90; -90, 90; -90, 90];

fprintf("Itération %d:\n", nb_iter_total);
fprintf("Masses initiales: [%f; %f; %f]\n", m_e_init);
fprintf("Masse totale initiale: %f\n", M);

% simulation de la trajectoire initiale
[R_tf, V_tf, y_1, y_2, y_3, t_1, t_2, t_3] = simulation_trajectoire(THETA_init(1), ...
            [THETA_init(2); THETA_init(3); THETA_init(4)], m_e_init, m_s, M, Mi, mu, c_x, R_terre, H_c, R_c, v_e, alpha);
        
fprintf("Vitesse réelle initiale = %f\n", norm(V_tf));

while 1
        nb_iter_total = nb_iter_total + 1;
        fprintf("\nItération %d:\n", nb_iter_total);

    %-----------------------------ETAGEMENT-------------------------------%
        V_p = V_p + deltaV;
    
        % Handle pour que ariane soit une fonction d'une variable
        ariane6 = @(x) ariane(x, V_p, m_satellite, k, v_e);
        fprintf("Vitesse propulsive = %f\n", V_p);       
        
        % Appel à SQP 
        [x_all, f_all, c_all, lambda_all, grad_L_norm_all, nb_iter, nb_eval, rho_all] = SQP(m_e_init, ...
        lambda_m_e_init, ariane6, @merite, 'SR1', bornes_m_e, 1e4, 1e7, 2, 1e-4, 2^12);
        fprintf("Nb d'itérations pb d'étagement : %d\n", nb_iter);
        
        % Masses optimales
        m_e_opt = x_all(:,end);
        fprintf("masses optimales: [%f; %f; %f]\n", m_e_opt);
        
        % masse sèche 
        m_s = k.*m_e_opt;
        % M est la masse initiale du lanceur et Mi est un vecteur de la
        % masse du lanceur à l'allumage de l'étage 1,2 et 3 (kg) 
        [M, ~, ~, Mi] = ariane6(m_e_opt);
        
        fprintf("Masse totale: %f\n", M);
        
    %-----------------------------TRAJECTOIRE-----------------------------%
    %----------------------------(branchement)----------------------------%
        
        % handles pour que la simulation soit une fonction de THETA seulement
        % et de même pour le branchement
        simulation_handle = @(THETA) simulation_trajectoire(THETA(1), [THETA(2); THETA(3); THETA(4)], m_e_opt, m_s, M, Mi, mu, c_x, R_terre, H_c, R_c, v_e, alpha);
        branchement_handle = @(THETA) branchement(THETA, simulation_handle, R_c);
        
        % Appel à SQP
        [THETA_all, f_THETA_all, c_THETA_all, lambda_THETA_all, grad_L_norm_THETA_all, nb_iter_THETA, nb_eval_THETA] = SQP(THETA_init, lambda_THETA_init, branchement_handle, @merite, 'BFGS', bornes_THETA, 1e3, 1e7, 2, 1e-8, 2^12);
        fprintf("Nb d'itérations pb de trajectoire : %d\n", nb_iter_THETA);
  
        % Theta optimale
        theta_opt = THETA_all(:,end);
        fprintf("theta optimal = [%f; %f; %f; %f]\n", theta_opt);
        
        % simulation de la trajectoire optimale
        [R_tf, V_tf, y_1, y_2, y_3, t_1, t_2, t_3] = simulation_trajectoire(theta_opt(1), ...
            [theta_opt(2); theta_opt(3); theta_opt(4)], m_e_opt, m_s, M, Mi, mu, c_x, R_terre, H_c, R_c, v_e, alpha);
        
        fprintf("Vitesse réelle = %f\n", norm(V_tf));

    % mise à jour de deltaV    
    deltaV = V_c - norm(V_tf);
    
    fprintf("deltaV = %f\n", deltaV);
    fprintf("norm(R_tf) - R_c = %f\n", norm(R_tf) - R_c);
    fprintf("R_tf'*V_tf = %f\n", R_tf'*V_tf);
    THETA_init = theta_opt;
    %m_e_init = m_e_opt;
    
    if (abs(deltaV) < 1e-3 && (abs(R_tf'*V_tf) < 1e-3))
        break;
    end
end

traces_simulation_trajectoire(y_1, y_2, y_3, t_1, t_2, t_3, R_terre, H_c, R_c, V_c);
fprintf("Nombre d'itérations au total : %d\n", nb_iter_total);
