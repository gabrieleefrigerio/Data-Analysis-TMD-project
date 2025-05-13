function X_hat = SystemSolver(i_nat, omega, density, A, E, J, M_a, L)
% SYSTEMSOLVER Calcola i coefficienti modali (X_hat) per ciascuna frequenza naturale.
% Risolve il sistema di equazioni associato alle condizioni al contorno
% per ricavare la forma generale della soluzione modale della trave.
%
% INPUT:
%   i_nat   - Indici delle frequenze naturali nel vettore omega
%   omega   - Vettore delle pulsazioni candidate
%   density - Densità del materiale della trave (ρ)
%   A       - Area della sezione trasversale della trave
%   E       - Modulo di Young del materiale
%   J       - Momento d'inerzia della sezione
%   M_a     - Massa aggiunta in estremità
%   L       - Lunghezza della trave
%
% OUTPUT:
%   X_hat   - Matrice 4xN dei coefficienti modali per ogni modo

    % Prealloca la matrice dei coefficienti modali: 4 righe (termini A, B, C, D), una colonna per ciascun modo
    X_hat = zeros(4, length(i_nat));
    
    % Cicla su ciascun modo proprio individuato
    for i_mode = 1:length(i_nat)
        % Estrae la pulsazione naturale corrente
        w = omega(i_nat(i_mode));
        
        % Calcola i parametri dinamici del sistema
        gamma = sqrt(w) * (density * A / (E * J))^(1/4);     % Parametro gamma
        tau = E * J * gamma^3;                                % Coefficiente legato alla rigidezza e inerzia
        lambda = w^2 * M_a;                                   % Effetto della massa concentrata
        gL = gamma * L;                                       % Gamma moltiplicato per la lunghezza (argomenti delle funzioni)

        % Costruisce la matrice del sistema (condizioni al contorno)
        Ab = [
            0, 1, 0, 1;                                                       % Condizione: w(0) = 0
            1, 0, 1, 0;                                                       % Condizione: w'(0) = 0
            -sin(gL), -cos(gL), sinh(gL), cosh(gL);                           % Condizione: momento flettente in L
            -tau*cos(gL)+lambda*sin(gL), ...                                  % Condizione: equilibrio forze a L
             tau*sin(gL)+lambda*cos(gL), ...
             tau*cosh(gL)+lambda*sinh(gL), ...
             tau*sinh(gL)+lambda*cosh(gL)
        ];

        % Isola il sistema lineare per risolvere i coefficienti B, C, D assumendo A = 1
        Hi_hat = Ab(2:4, 2:4);  % Matrice dei coefficienti sconosciuti
        Ni_hat = Ab(2:4, 1);    % Parte nota del sistema (prima colonna con A=1)

        % Risolve il sistema lineare: X_hat = [A; B; C; D] con A = 1 fissato
        X_hat(:, i_mode) = [1; -Hi_hat \ Ni_hat];
    end
end

