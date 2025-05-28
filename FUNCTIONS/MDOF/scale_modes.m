function scaled_modes = scale_modes(FRF, modal_results)
% scale_modes - Scala le forme modali tramite least-squares sulla FRF

    freq = modal_results.freq;
    modes = modal_results.modes;
    poles = modal_results.poles;
    % Converti freq [Hz] → omega [rad/s]
    omega = 2 * pi * FRF.freq;  % vettore colonna
    Nf = length(omega);
    Nm = length(modes);  % numero di modi

    % Assumiamo H_{pq} = FRF(:, p, q), seleziona p=2, q=1 (modifica se necessario)
    H_measured = squeeze(FRF.FRF(:, 2, 1));  % vettore colonna complesso

    % Prepara matrice A: Nf righe, Nm+2 colonne
    A = zeros(Nf, Nm + 2);

    for r = 1:Nm
        psi = modes{r};  % vettore colonna
        if isempty(psi)
            continue;
        end

        % Supponiamo il prodotto ψ_{pr} * ψ_{qr} come primo elemento del modo (semplificazione)
        % Altrimenti definisci p, q esplicitamente
        psi_prod = psi(1) * psi(2);  % esempio p=1, q=2

        s_r = poles(r);  % polo complesso
        A(:, r) = psi_prod ./ (1i * omega - s_r);  % contributo modale
    end

    % Colonne residue: R_U (constante), R_L (1/omega^2)
    A(:, Nm+1) = 1;              % R_U
    A(:, Nm+2) = 1 ./ omega.^2;  % R_L

    % Risolvi sistema: A x ≈ H_measured
    x = A \ H_measured;

    % Estrai Q_r
    Q = real(x(1:Nm));  % forza residua associata a ciascun modo

    % Scala le forme modali: φ_scaled = ψ * sqrt(Q_r)
    scaled_modes = cell(Nm, 1);
    for r = 1:Nm
        if Q(r) > 0 && ~isempty(modes{r})
            scaled_modes{r} = modes{r} * sqrt(Q(r));
        else
            scaled_modes{r} = [];
        end
    end
end
