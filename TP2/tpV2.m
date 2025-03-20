close all;
clear all;
clc;

% Création du dossier pour stocker les figures
output_folder = 'figures';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Chargement des matrices disponibles
matrices = {'mat0', 'mat1', 'mat2', 'mat3'};
permutation_methods = {'RCM', 'AMD'};
resultats = [];

for i = 1:length(matrices)
    load(matrices{i}); % Charger la matrice
    fprintf('Traitement de la matrice : %s\n', matrices{i});

    n = size(A,1);
    b = (1:n)';

    % Boucle sur les différentes méthodes de permutation
    for j = 1:length(permutation_methods)
        method = permutation_methods{j};
        
        % Sélection de la permutation
        switch method
            case 'RCM'
                P = symrcm(A); % Reverse Cuthill-McKee
            case 'AMD'
                P = amd(A); % Approximate Minimum Degree
        end

        % Matrice permutée
        B = A(P,P);

        % Création d'une figure unique pour chaque méthode et matrice
        figure;
        
        % Matrice originale
        subplot(2,3,1);
        spy(A);
        title('Original matrix A');

        [count, h, parent, post, R] = symbfact(A);
        ALU = R + R';
        subplot(2,3,2);
        spy(ALU);
        title('Factors of A');
        fillin_no_perm = nnz(ALU) - nnz(A);

        % Visualisation du remplissage
        C = spones(A);
        CLU = spones(ALU);
        FILL = CLU - C;
        subplot(2,3,3);
        spy(FILL);
        title('Fill on original A');

        % Matrice permutée
        subplot(2,3,4);
        spy(B);
        title(['Permuted matrix (', method, ')']);

        [count, h, parent, post, R] = symbfact(B);
        BLU = R + R';
        subplot(2,3,5);
        spy(BLU);
        fillin_perm = nnz(BLU) - nnz(A);
        title(['Factors of permuted A (', method, ')']);

        % Visualisation du remplissage après permutation
        B = spones(B);
        BLU = spones(BLU);
        FILL = BLU - B;
        subplot(2,3,6);
        spy(FILL);
        title(['Fill on permuted A (', method, ')']);

        % Sauvegarde de la figure
        saveas(gcf, fullfile(output_folder, sprintf('%s_%s.png', matrices{i}, method)));

        % Résolution sans réordonnancement
        x_ref = A \ b;
        L = chol(A, 'lower');
        flops_no_perm = 4 * nnz(L);
        y = L \ b;
        x = L' \ y;
        error_no_perm = norm(x - x_ref) / norm(x_ref);

        % Résolution avec réordonnancement
        B = A(P,P);
        L = chol(B, 'lower');
        flops_perm = 4 * nnz(L);
        b_perm = b(P);
        y = L \ b_perm;
        x_perm = L' \ y;
        x(P) = x_perm;
        error_perm = norm(x - x_ref) / norm(x_ref);

        % Stockage des résultats
        resultats = [resultats; {matrices{i}, method, fillin_no_perm, fillin_perm, flops_no_perm, flops_perm, error_no_perm, error_perm}];

        % Fermer la figure
        close;
    end
end

% Sauvegarde des résultats
T = cell2table(resultats, 'VariableNames', {'Matrice', 'Methode', 'Fill-in_Sans', 'Fill-in_Avec', 'Flops_Sans', 'Flops_Avec', 'Erreur_Sans', 'Erreur_Avec'});
writetable(T, 'resultats_comparaison.csv');
