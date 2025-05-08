function G = TransferFunctionComputation(data)
    %TRANSFERFUNCTIONCOMPUTATION Calcola le funzioni di trasferimento (FRF),
    % le frequenze naturali e le forme modali di una trave cilindrica cava
    % a sbalzo con massa concentrata in punta, tramite metodo modale.
    
    %% === Parametri geometrici e materiali ===
    rho   = data.rho;             % Densità materiale [kg/m^3]
    L     = data.L;               % Lunghezza della trave [m]
    E     = data.Ea * 1e9;        % Modulo di Young [Pa]
    D_ext = data.D_ext;           % Diametro esterno [m]
    D_int = data.D_int;           % Diametro interno [m]
    M_a   = data.M;               % Massa concentrata in punta [kg]
    
    % Area e momento d'inerzia della sezione trasversale
    A = (pi / 4) * (D_ext^2 - D_int^2);         
    J = (pi / 64) * (D_ext^4 - D_int^4);        
    
    %% === Frequenze di calcolo ===
    fmax = 80;                              % Frequenza massima [Hz]
    n_points = 10000;                       % Numero di punti in frequenza
    freq = linspace(0, fmax, n_points);    % Vettore frequenze [Hz]
    omega = 2 * pi * freq;                 % Frequenze angolari [rad/s]
    
    %% === Calcolo del determinante della matrice e frequenze naturali ===
    dets = MatrixDeterminant(omega, rho, A, E, J, M_a, L);        
    [i_nat, G] = NaturalFrequencyComputation(dets, freq);       
    
    %% === Calcolo coefficienti modali ===
    X_hat = SystemSolver(i_nat, omega, rho, A, E, J, M_a, L);
    
    %% === Costruzione delle forme modali ===
    [x, modes_shapes] = ModeShapes(L, n_points, i_nat, omega, rho, A, E, J, X_hat);
    
    %% === Definizione posizione forza e sensore ===
    xj = data.L;  % Posizione accelerometro [m]
    xk = data.L;  % Posizione applicazione forza [m]
    [~, pos_xj] = min(abs(x - xj));  
    [~, pos_xk] = min(abs(x - xk));  % (non usato direttamente)
    
    %% === Calcolo massa modale ===
    m = rho * A;  
    modal_mass = zeros(1, length(i_nat));
    for i = 1:length(i_nat)
        modal_mass(i) = m * trapz(x, modes_shapes(i,:).^2);
    end
    
    %% === Calcolo delle FRF: G_WY e G_SY ===
    xsi = data.xi1;  % Smorzamento critico (es: 0.01)
    [phi_at_xj, G] = GwyComputation(omega, i_nat, modes_shapes, x, pos_xj, modal_mass, xsi, A, M_a, G);
    G = GsyComputation(omega, i_nat, modal_mass, xsi, phi_at_xj, G);
    
    %% === Plot finale ===
    GPlotSingle(freq, G);  % Visualizzazione risultati
end

%% AUXILIAR FUNCTIONS

function dets = MatrixDeterminant(omega, rho, A, E, J, M_a, L)
% Calcola il determinante della matrice del sistema per ciascuna frequenza.
    
    dets = zeros(length(omega), 1);
    for i = 1:length(omega)
        w = omega(i);
        gamma = sqrt(w) * (rho * A / (E * J))^(1/4);
        tau = E * J * gamma^3;
        lambda = w^2 * M_a;
        gL = gamma * L;
    
        Ab = [
            0, 1, 0, 1;
            1, 0, 1, 0;
            -sin(gL), -cos(gL), sinh(gL), cosh(gL);
            -tau * cos(gL) + lambda * sin(gL), ...
             tau * sin(gL) + lambda * cos(gL), ...
             tau * cosh(gL) + lambda * sinh(gL), ...
             tau * sinh(gL) + lambda * cosh(gL)
        ];
    
        dets(i) = det(Ab);
    end
end

function [i_nat, G] = NaturalFrequencyComputation(dets, freq)
% Trova gli indici delle frequenze naturali dai minimi locali del determinante.

    i_nat = [];
    for i = 2:length(dets)-1
        if abs(dets(i)) < abs(dets(i-1)) && abs(dets(i)) < abs(dets(i+1))
            i_nat(end+1) = i;
        end
    end
    G.freq_nat = freq(i_nat);
    end
        
function X_hat = SystemSolver(i_nat, omega, rho, A, E, J, M_a, L)
    % Risolve il sistema per ottenere i coefficienti modali X_hat per ogni modo.
    
    X_hat = zeros(4, length(i_nat));
    for i_mode = 1:length(i_nat)
        w = omega(i_nat(i_mode));
        gamma = sqrt(w) * (rho * A / (E * J))^(1/4);
        tau = E * J * gamma^3;
        lambda = w^2 * M_a;
        gL = gamma * L;
    
        Ab = [
            0, 1, 0, 1;
            1, 0, 1, 0;
            -sin(gL), -cos(gL), sinh(gL), cosh(gL);
            -tau * cos(gL) + lambda * sin(gL), ...
             tau * sin(gL) + lambda * cos(gL), ...
             tau * cosh(gL) + lambda * sinh(gL), ...
             tau * sinh(gL) + lambda * cosh(gL)
        ];
    
        Hi_hat = Ab(2:4, 2:4);
        Ni_hat = Ab(2:4, 1);
        X_hat(:, i_mode) = [1; -Hi_hat \ Ni_hat];
    end
end

function [x, modes_shapes] = ModeShapes(L, n_points, i_nat, omega, rho, A, E, J, X_hat)
% Costruisce le forme modali normalizzate lungo la trave.

    x = linspace(0, L, n_points);
    modes_shapes = zeros(length(i_nat), length(x));
    
    for i_mode = 1:length(i_nat)
        w = omega(i_nat(i_mode));
        gamma = sqrt(w) * (rho * A / (E * J))^(1/4);
    
        modes_shapes(i_mode, :) = ...
            X_hat(1, i_mode) * sin(gamma * x) + ...
            X_hat(2, i_mode) * cos(gamma * x) + ...
            X_hat(3, i_mode) * sinh(gamma * x) + ...
            X_hat(4, i_mode) * cosh(gamma * x);
    end
    
    % Normalizzazione rispetto al valore massimo
    modes_shapes = modes_shapes ./ max(abs(modes_shapes), [], 2);
end

function [phi_at_xj, G] = GwyComputation(omega, i_nat, modes_shapes, x, pos_xj, modal_mass, xsi, A, M_a, G)
% Calcola la funzione di trasferimento teorica G_WY (accelerazione).

    G_FRF_theory = zeros(length(omega), 1);
    phi_integrals = zeros(1, length(i_nat));
    phi_at_xj = zeros(1, length(i_nat));
    phiL = zeros(1, length(i_nat));
    
    for i = 1:length(i_nat)
        phi_i = modes_shapes(i, :);
        phi_integrals(i) = trapz(x, phi_i);
        phi_at_xj(i) = phi_i(pos_xj);
        phiL(i) = phi_i(end);
    end
    
    for k = 1:length(omega)
        w = omega(k);
        sum_modes = 0;
        for i = 1:length(i_nat)
            w_n = omega(i_nat(i));
            mi = modal_mass(i);
            zeta_i = xsi;
            numerator = w^2 * phi_at_xj(i);
            denominator = mi * (-w^2 + 2j * zeta_i * w * w_n + w_n^2);
            participation = A * phi_integrals(i) + M_a * phiL(i);
            sum_modes = sum_modes + (numerator / denominator) * participation;
        end
        G_FRF_theory(k) = 1 + sum_modes;
    end
    
    G.Gwy = G_FRF_theory;
    G.Gwy_amp = abs(G_FRF_theory);
    G.Gwy_phase = angle(G_FRF_theory);
end

function G = GsyComputation(omega, i_nat, modal_mass, xsi, phi_at_xj, G)
% Calcola la funzione di trasferimento numerica G_SY (spostamento).

    G_sy = zeros(length(omega), 1);
    for k = 1:length(omega)
        w = omega(k);
        sum_modes_sy = 0;
        for i = 1:length(i_nat)
            w_n = omega(i_nat(i));
            mi = modal_mass(i);
            zeta_i = xsi;
            numerator = phi_at_xj(i);
            denominator = mi * (-w^2 + 2j * zeta_i * w * w_n + w_n^2);
            sum_modes_sy = sum_modes_sy + (numerator / denominator);
        end
        G_sy(k) = sum_modes_sy;
    end
    
    G.Gsy = G_sy;
    G.Gsy_amp = abs(G_sy);
    G.Gsy_phase = angle(G_sy);
end

function GPlotSingle(freq, G)
% Crea i grafici separati di ampiezza e fase per G_WY e G_SY

    % G_WY - Ampiezza
    figure('Color','w','WindowStyle','docked');
    semilogy(freq, abs(G.Gwy), 'b', 'LineWidth', 1.6);
    grid on; xlabel('Frequenza [Hz]'); ylabel('|G_{WY}|'); title('Ampiezza G_{WY}');
    
    % G_WY - Fase
    figure('Color','w','WindowStyle','docked');
    plot(freq, unwrap(angle(G.Gwy))*180/pi, 'b', 'LineWidth', 1.6);
    grid on; xlabel('Frequenza [Hz]'); ylabel('Fase G_{WY} [°]'); title('Fase G_{WY}');
    
    % G_SY - Ampiezza
    figure('Color','w','WindowStyle','docked');
    semilogy(freq, abs(G.Gsy), 'r', 'LineWidth', 1.6);
    grid on; xlabel('Frequenza [Hz]'); ylabel('|G_{SY}|'); title('Ampiezza G_{SY}');
    
    % G_SY - Fase
    figure('Color','w','WindowStyle','docked');
    plot(freq, unwrap(angle(G.Gsy))*180/pi, 'r', 'LineWidth', 1.6);
    grid on; xlabel('Frequenza [Hz]'); ylabel('Fase G_{SY} [°]'); title('Fase G_{SY}');
end

