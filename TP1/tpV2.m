close all;
clear all;


%load mat1;
%load pde225_5e-1;
load hydcar20.mat;
n = size(A, 1);
fprintf('Dimension de A : %4d \n', n);


b = [1:n]';   % Le vecteur b de taille n
x0 = zeros(n, 1); % Point de départ initial
kmax = n-1;     % Nombre maximum d'itérations
A_dense = full(A);

% Calculer le nombre de conditionnement de A_dense
condition_number = cond(A_dense);
fprintf('Nombre de conditionnement de A : %f\n', condition_number);

% Initialisation d'un tableau pour stocker les résultats
results = table(); % Table vide pour stocker les résultats
iteration_count = 0; % Compteur d'itérations pour les résultats

%% Test pour différentes tolérances
tolerances = [1e-6, 1e-8, 1e-10];

for tol_idx = 1:length(tolerances)
    eps = tolerances(tol_idx);

    % FOM
    [x_fom, flag_fom, relres_fom, iter_fom, resvec_fom] = krylov(A, b, x0, eps, kmax, 0);
    iteration_count = iteration_count + 1;
    results(iteration_count, :) = {eps, 'FOM', iter_fom, relres_fom, flag_fom};

    % GMRES
    [x_gmres, flag_gmres, relres_gmres, iter_gmres, resvec_gmres] = krylov(A, b, x0, eps, kmax, 1);
    iteration_count = iteration_count + 1;
    results(iteration_count, :) = {eps, 'GMRES', iter_gmres, relres_gmres, flag_gmres};

    % GMRES MATLAB
    [x_gmres_matlab, flag_gmres_matlab, relres_gmres_matlab, iter_gmres_matlab, resvec_gmres_matlab] = gmres(A, b, [], eps, kmax, [], [], x0);
    iteration_count = iteration_count + 1;
    results(iteration_count, :) = {eps, 'GMRES MATLAB', iter_gmres_matlab(2), relres_gmres_matlab, flag_gmres_matlab};

    % Graphiques
    figure;
    semilogy(resvec_fom, 'c', 'DisplayName', 'FOM');
    hold on;
    semilogy(resvec_gmres, 'r', 'DisplayName', 'GMRES');
    semilogy(resvec_gmres_matlab, 'b+', 'DisplayName', 'GMRES MATLAB');
    legend('show');
    title(['Évolution de la norme du résidu pour \epsilon = ', num2str(eps)]);
    xlabel('Nombre d\ itérations');
    ylabel('Norme du résidu');

    % Enregistrement de la figure
    saveas(gcf, ['residu_epsilon_' num2str(eps) '.png']);
end

% Enregistrement des résultats
writetable(results, 'krylov_results.csv');
