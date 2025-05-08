function G = TransferFunctionComputation(data)
% TransferFunctionComputation
% Calcola frequenze naturali e modi propri di una trave cilindrica cava in cantilever
% con massa concentrata in punta, usando il metodo della funzione di trasferimento.

    %% --- Parametri geometrici e materiali ---
    rho = data.rho;             % densità [kg/m^3]
    L = data.L;                 % lunghezza trave [m]
    E = data.Ea * 1e9;          % modulo di Young [Pa]
    D_ext = data.D_ext;        % diametro esterno [m]
    D_int = data.D_int;        % diametro interno [m]
    M_a = data.M;              % massa concentrata in punta [kg]

    A = (pi / 4) * (D_ext^2 - D_int^2);        % area [m^2]
    J = (pi / 64) * (D_ext^4 - D_int^4);       % momento d'inerzia [m^4]
    
    %% --- Asse delle frequenze ---
    fmax = 200;               % max frequenza [Hz]
    n_points = 10000;
    freq = linspace(0, fmax, n_points);
    omega = 2 * pi * freq;

    %% --- Calcolo determinante della matrice Ab ---
    dets = zeros(length(omega), 1);

    for i = 1:length(omega)
        w = omega(i);
        gamma = sqrt(w) * (rho * A / (E * J))^(1/4);   % formula modificata
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

    %% --- Ricerca delle frequenze naturali ---
    i_nat = [];
    for i = 2:length(dets)-1
        if abs(dets(i)) < abs(dets(i-1)) && abs(dets(i)) < abs(dets(i+1))
            i_nat(end+1) = i;
        end
    end
    G.freq_nat = freq(i_nat);

    %% --- Calcolo coefficienti modali X_hat ---
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

    %% --- Costruzione forme modali ---
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

    % Normalizzazione
    modes_shapes = modes_shapes ./ max(abs(modes_shapes), [], 2);

    %% --- Output struttura ---
    G.freq = freq;
    G.detH = dets;
    G.X_hat = X_hat;
    G.modes_shapes = modes_shapes;

    %% --- (Facoltativa) Animazione forme modali ---
    
    mode = 1;  % seleziona il modo
    T = 1 / G.freq_nat(mode);
    n_frames = 500;
    t = linspace(0, 4*T, n_frames);

    figure('Position', [100, 100, 1200, 400]); hold on;
    title(['Animazione - Modo ', num2str(mode)])
    plot(x, G.modes_shapes(mode,:), '--k', 'LineWidth', 2)
    h1 = plot(x, zeros(size(x)), 'LineWidth', 2);
    xlabel('Posizione [m]'); ylabel('Ampiezza normalizzata')
    ylim([-1.1 1.1]);

    for k = 1:length(t)
        h1.YData = G.modes_shapes(mode,:) * cos(2*pi*G.freq_nat(mode)*t(k));
        pause(T / n_frames);
    end
    




    % %% Transfer Function
    % 
    % % posizioni acceleroemtri
    % xj = [ 0.1 0.2 0.3 0.5 0.7 0.8 1 1.2]; %[m]
    % 
    % % posizione martellata
    % xk = 1.2; %[m]
    % 
    % 
    % % Inizializzazione matrice FRF per tutti gli xj
    % frf = zeros(length(omega), length(xj));
    % 
    % % Colori per i plot (opzionale: si può personalizzare con colormap)
    % color_map = lines(length(xj));
    % 
    % % Calcolo della FRF per ogni xj
    % for ii = 1:length(xj)
    %     % Trova gli indici corrispondenti a xj e xk
    %     [~, pos_xj] = min(abs(x - xj(ii)));
    %     [~, pos_xk] = min(abs(x - xk));
    % 
    %     % Calcolo massa modale
    %     modal_mass = zeros(1, length(i_nat));
    %     for i = 1:length(i_nat)
    %         modal_mass(i) = m * trapz(x, modes_shapes(i,:).^2);
    %     end
    % 
    %     % Calcolo FRF
    %     FRF = zeros(1, length(omega));
    %     for k = 1:length(omega)
    %         for i = 1:length(i_nat)
    %             num = modes_shapes(i,pos_xj) * modes_shapes(i,pos_xk);
    %             den = -omega(k)^2 + 2i*xsi*omega(k)*omega(i_nat(i)) + omega(i_nat(i))^2;
    %             FRF(k) = FRF(k) + (num / modal_mass(i)) / den;
    %         end
    %     end
    % 
    %     % Salva FRF calcolata
    %     frf(:, ii) = FRF.';
    % end
    % 
    % % === PLOT ===
    % 
    % % Calcolo ampiezza e fase
    % FRF_amp = abs(frf);
    % FRF_phase = angle(frf);
    % 
    % % Crea la figura
    % figure('Name', 'FRF - Ampiezza e Fase', 'Color', 'w');
    % 
    % % ---- Subplot Ampiezza (semilog) ----
    % subplot(2,1,1);
    % hold on;
    % for ii = 1:length(xj)
    %     semilogy(freq, FRF_amp(:,ii), 'LineWidth', 1.8, 'Color', color_map(ii,:));
    %     set(gca, 'YScale', 'log')
    % end
    % grid on;
    % xlabel('Frequenza [Hz]', 'FontSize', 11);
    % ylabel('|G(j\omega)|', 'FontSize', 11);
    % title(' Ampiezza', 'FontWeight', 'bold');
    % xlim([min(freq) max(freq)]);
    % legend(arrayfun(@(v) sprintf('acc in %.1f m', v), xj, 'UniformOutput', false), ...
    %        'FontSize', 9, 'Location', 'northeast');
    % 
    % % ---- Subplot Fase ----
    % subplot(2,1,2);
    % hold on;
    % for ii = 1:length(xj)
    %     plot(freq, FRF_phase(:,ii), 'LineWidth', 1.8, 'Color', color_map(ii,:));
    % end
    % grid on;
    % xlabel('Frequenza [Hz]', 'FontSize', 11);
    % ylabel('Fase [rad]', 'FontSize', 11);
    % title('Fase', 'FontWeight', 'bold');
    % xlim([min(freq) max(freq)]);
    % legend(arrayfun(@(v) sprintf('acc in %.1f m', v), xj, 'UniformOutput', false), ...
    %        'FontSize', 9, 'Location', 'northeast');
    % 
    

end

