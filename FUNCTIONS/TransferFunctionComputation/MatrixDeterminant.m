function dets = MatrixDeterminant(omega, density, A, E, J, M_a, L)
% MATRIXDETERMINANT Calcola il determinante della matrice del sistema H(ω)
% per ciascuna frequenza naturale fornita.
%
% INPUT:
%   omega   - Vettore delle frequenze angolari [rad/s]
%   density - Densità del materiale [kg/m^3]
%   A       - Area della sezione trasversale [m^2]
%   E       - Modulo di Young [Pa]
%   J       - Momento d'inerzia della sezione [m^4]
%   M_a     - Massa concentrata in estremità [kg]
%   L       - Lunghezza della trave [m]
%
% OUTPUT:
%   dets    - Vettore dei determinanti calcolati per ogni ω

    % Preallocazione vettore output
    dets = zeros(length(omega), 1);

    % Loop su ogni frequenza ω
    for i = 1:length(omega)
        w = omega(i);  % Frequenza corrente

        % Calcolo parametri dipendenti dalla frequenza
        gamma = sqrt(w) * (density * A / (E * J))^(1/4);  % Parametro gamma
        tau = E * J * gamma^3;                            % Costante tau
        lambda = w^2 * M_a;                               % Parametro lambda
        gL = gamma * L;                                   % Prodotto gamma * L

        % Matrice H(w): matrice del sistema per il problema di vincolo
        Hmat = [
            0,                          1,                          0,                          1;
            1,                          0,                          1,                          0;
           -sin(gL),                  -cos(gL),                   sinh(gL),                   cosh(gL);
           -tau*cos(gL) + lambda*sin(gL),  tau*sin(gL) + lambda*cos(gL), ...
            tau*cosh(gL) + lambda*sinh(gL), tau*sinh(gL) + lambda*cosh(gL)
        ];

        % Calcolo del determinante della matrice
        dets(i) = det(Hmat);
    end
end


