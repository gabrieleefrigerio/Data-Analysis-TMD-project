function G = TransferFunctionComputation(data)
    % TransferFunctionComputation
    % Calcola frequenze naturali, modi propri e funzioni di trasferimento
    % per una trave cilindrica cava in condizioni di vincolo a sbalzo
    % con massa concentrata in punta, utilizzando il metodo della funzione di trasferimento.

    %% --- Parametri geometrici e materiali ---
    % Definizione dei parametri geometrici e dei materiali
    rho = data.rho;             % densità [kg/m^3]
    L = data.L;                 % lunghezza trave [m]
    E = data.Ea * 1e9;          % modulo di Young [Pa]
    D_ext = data.D_ext;         % diametro esterno [m]
    D_int = data.D_int;         % diametro interno [m]
    M_a = data.M;               % massa concentrata in punta [kg]

    % Calcoli dell'area e del momento di inerzia
    A = (pi / 4) * (D_ext^2 - D_int^2);         % area [m^2]
    J = (pi / 64) * (D_ext^4 - D_int^4);        % momento d'inerzia [m^4]
    
    %% --- Asse delle frequenze ---
    % Frequenze di interesse per il calcolo della funzione di trasferimento
    fmax = 80;              % Frequenza massima [Hz]
    n_points = 10000;       % Numero di punti sulla griglia delle frequenze
    freq = linspace(0, fmax, n_points);  % Vettore delle frequenze [Hz]
    omega = 2 * pi * freq;  % Frequenze angolari [rad/s]

    %% --- Calcolo determinante della matrice Ab ---
    % Determinante della matrice Ab per trovare le frequenze naturali
    dets = zeros(length(omega), 1);  % Vettore per i determinanti

    for i = 1:length(omega)
        w = omega(i);  % Frequenza angolare
        gamma = sqrt(w) * (rho * A / (E * J))^(1/4);  % Parametro gamma
        tau = E * J * gamma^3;
        lambda = w^2 * M_a;
        gL = gamma * L;

        % Costruzione della matrice Ab
        Ab = [
            0, 1, 0, 1;
            1, 0, 1, 0;
            -sin(gL), -cos(gL), sinh(gL), cosh(gL);
            -tau * cos(gL) + lambda * sin(gL), ...
             tau * sin(gL) + lambda * cos(gL), ...
             tau * cosh(gL) + lambda * sinh(gL), ...
             tau * sinh(gL) + lambda * cosh(gL)
        ];

        % Calcolo del determinante
        dets(i) = det(Ab);
    end

    %% --- Ricerca delle frequenze naturali ---
    % Identifica le frequenze naturali da dove il determinante cambia segno
    i_nat = [];
    for i = 2:length(dets)-1
        if abs(dets(i)) < abs(dets(i-1)) && abs(dets(i)) < abs(dets(i+1))
            i_nat(end+1) = i;  % Aggiungi l'indice della frequenza naturale
        end
    end
    G.freq_nat = freq(i_nat);  % Frequenze naturali

    %% --- Calcolo dei coefficienti modali X_hat ---
    % Calcolo dei coefficienti modali per ogni modo naturale
    X_hat = zeros(4, length(i_nat));

    for i_mode = 1:length(i_nat)
        w = omega(i_nat(i_mode));
        gamma = sqrt(w) * (rho * A / (E * J))^(1/4);
        tau = E * J * gamma^3;
        lambda = w^2 * M_a;
        gL = gamma * L;

        % Costruzione della matrice Ab
        Ab = [
            0, 1, 0, 1;
            1, 0, 1, 0;
            -sin(gL), -cos(gL), sinh(gL), cosh(gL);
            -tau * cos(gL) + lambda * sin(gL), ...
             tau * sin(gL) + lambda * cos(gL), ...
             tau * cosh(gL) + lambda * sinh(gL), ...
             tau * sinh(gL) + lambda * cosh(gL)
        ];

        % Calcolo dei coefficienti modali X_hat
        Hi_hat = Ab(2:4, 2:4);
        Ni_hat = Ab(2:4, 1);
        X_hat(:, i_mode) = [1; -Hi_hat \ Ni_hat];
    end

    %% --- Costruzione delle forme modali ---
    % Costruzione delle forme modali della trave
    x = linspace(0, L, n_points);
    modes_shapes = zeros(length(i_nat), length(x));

    for i_mode = 1:length(i_nat)
        w = omega(i_nat(i_mode));
        gamma = sqrt(w) * (rho * A / (E * J))^(1/4);

        % Calcolo delle forme modali
        modes_shapes(i_mode, :) = ...
            X_hat(1, i_mode) * sin(gamma * x) + ...
            X_hat(2, i_mode) * cos(gamma * x) + ...
            X_hat(3, i_mode) * sinh(gamma * x) + ...
            X_hat(4, i_mode) * cosh(gamma * x);
    end

    % Normalizzazione delle forme modali
    modes_shapes = modes_shapes ./ max(abs(modes_shapes), [], 2);

    %% --- Output della struttura ---
    % Restituisce i risultati calcolati
    G.freq = freq;              % Frequenze [Hz]
    G.detH = dets;              % Determinanti della matrice Ab
    G.X_hat = X_hat;            % Coefficienti modali
    G.modes_shapes = modes_shapes;  % Forme modali

    %% --- Calcolo della Funzione di Trasferimento (FRF) ---
    % Definizione delle posizioni degli accelerometri e della forza
    xj = data.L;  % Accelerometro in punta
    xk = data.L;  % Forza applicata in punta

    % Trova gli indici corrispondenti sulle griglie x
    [~, pos_xj] = min(abs(x - xj));
    [~, pos_xk] = min(abs(x - xk));

    % Calcolo della massa modale per ogni modo
    m = rho * A;  % Massa per unità di lunghezza
    modal_mass = zeros(1, length(i_nat));
    for i = 1:length(i_nat)
        modal_mass(i) = m * trapz(x, modes_shapes(i,:).^2);  % Massa modale
    end

    % Fattore di smorzamento modale
    xsi = data.xi1;  % Fattore di smorzamento critico (1%)

    %% --- Calcolo della FRF teorica G_WY ---
    % Calcola la funzione di trasferimento G_WY in base alla teoria
    theta_m = 1;  % Coefficiente di partecipazione modale
    G_FRF_theory = zeros(length(omega), 1);

    % Integrali delle forme modali
    phi_integrals = zeros(1, length(i_nat));
    phi_at_xj = zeros(1, length(i_nat));
    phiL = zeros(1, length(i_nat));

    for i = 1:length(i_nat)
        phi_i = modes_shapes(i, :);
        phi_integrals(i) = trapz(x, phi_i);    % ∫ϕᵢ(x)dx
        phi_at_xj(i) = phi_i(pos_xj);          % ϕᵢ(xj)
        phiL(i) = phi_i(end);                   % ϕᵢ(L)
    end

    % Calcolo della FRF teorica per ogni frequenza
    for k = 1:length(omega)
        w = omega(k);
        sum_modes = 0;

        for i = 1:length(i_nat)
            w_n = omega(i_nat(i));
            mi = modal_mass(i);
            zeta_i = xsi;

            % Calcolo del numeratore e denominatore per ogni modo
            numerator = w^2 * phi_at_xj(i);
            denominator = mi * (-w^2 + 2j * zeta_i * w * w_n + w_n^2);
            participation = theta_m * A * phi_integrals(i) + M_a * phiL(i);

            % Somma sui modi
            sum_modes = sum_modes + (numerator / denominator) * participation;
        end

        G_FRF_theory(k) = 1 + sum_modes;  % FRF teorica
    end

    % Risultati della FRF teorica
    G.Gwy = G_FRF_theory;
    G.Gwy_amp = abs(G_FRF_theory);
    G.Gwy_phase = angle(G_FRF_theory);

    %% --- Calcolo di Gsy ---
    % Calcolo della funzione di trasferimento G_sy (risposta di spostamento)
    G_sy = zeros(length(omega), 1);

    for k = 1:length(omega)
        w = omega(k);
        sum_modes_sy = 0;

        for i = 1:length(i_nat)
            w_n = omega(i_nat(i));
            mi = modal_mass(i);
            zeta_i = xsi;

            % Calcolo del numeratore e denominatore per ogni modo
            numerator = phi_at_xj(i);  % Risposta alla forza
            denominator = mi * (-w^2 + 2j * zeta_i * w * w_n + w_n^2);

            % Somma sui modi per G_sy
            sum_modes_sy = sum_modes_sy + (numerator / denominator);
        end

        G_sy(k) = sum_modes_sy;  % Funzione di trasferimento G_sy
    end

    % Risultati della FRF di spostamento
    G.Gsy = G_sy;
    G.Gsy_amp = abs(G_sy);
    G.Gsy_phase = angle(G_sy);

    %% --- Plot della FRF teorica vs numerica ---
    % Crea una figura separata per ogni grafico, così che possano essere dockabili
    % Colori personalizzati
    color_th = [0 0.4470 0.7410];  % Blu
    
    % Ampiezza della FRF - G_{WY}
    figure('Color', 'w', 'WindowStyle', 'docked');
    semilogy(freq, abs(G.Gwy), 'Color', color_th, 'LineWidth', 1.6);
    grid on;
    xlabel('Frequenza [Hz]');
    ylabel('|G_{WY}|');
    title('Ampiezza G_{WY}');
    legend('G_{WY}', 'Location', 'NorthEast');
    
    % Fase della FRF - G_{WY}
    figure('Color', 'w', 'WindowStyle', 'docked');
    plot(freq, unwrap(angle(G.Gwy)) * 180 / pi, 'Color', color_th, 'LineWidth', 1.6);
    grid on;
    xlabel('Frequenza [Hz]');
    ylabel('Fase G_{WY} [°]');
    title('Fase G_{WY}');
    legend('Fase G_{WY}', 'Location', 'Best');
    
    % Ampiezza della FRF - G_{SY}
    figure('Color', 'w', 'WindowStyle', 'docked');
    semilogy(freq, abs(G.Gsy), 'r', 'LineWidth', 1.6);
    grid on;
    xlabel('Frequenza [Hz]');
    ylabel('|G_{SY}|');
    title('Ampiezza G_{SY}');
    legend('G_{SY}', 'Location', 'NorthEast');
    
    % Fase della FRF - G_{SY}
    figure('Color', 'w', 'WindowStyle', 'docked');
    plot(freq, unwrap(angle(G.Gsy)) * 180 / pi, 'r', 'LineWidth', 1.6);
    grid on;
    xlabel('Frequenza [Hz]');
    ylabel('Fase G_{SY} [°]');
    title('Fase G_{SY}');
    legend('Fase G_{SY}', 'Location', 'Best');

end
