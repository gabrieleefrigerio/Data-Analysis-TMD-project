function G = GsyComputation(omega, i_nat, modes_shapes, modal_mass, xsi, phi_at_xj, G, E, J, x, pos_xj, density, A, M_a)
% GSYCOMPUTATION Calcola la funzione di trasferimento teorica G_SY (spostamento).

    G_sy = zeros(length(omega), 1);  % Inizializza il vettore di G_SY
    phi_integrals = zeros(1, length(i_nat));       
    phiL = zeros(1, length(i_nat));                
    d3_phi_at_0 = zeros(1, length(i_nat));         % Terza derivata in x=0

    dx = x(2) - x(1);  % Assumiamo passo uniforme

    % Calcolo dei coefficienti modali necessari (integrali, valori locali)
    for i = 1:length(i_nat)
        phi_i = modes_shapes(i, :);
        
        % Calcolo della derivata terza numerica in x=0
        % Usiamo differenze finite in avanti di ordine 3 (accuratezza O(h^2))
        % f'''(x0) ≈ ( -5f0 + 18f1 - 24f2 + 14f3 - 3f4 ) / (2h^3)
        if length(x) >= 5
            d3_phi_at_0(i) = (-5*phi_i(1) + 18*phi_i(2) - 24*phi_i(3) + 14*phi_i(4) - 3*phi_i(5)) / (2*dx^3);
        else
            error('Non ci sono abbastanza punti per calcolare la derivata terza a x=0');
        end

        phi_integrals(i) = trapz(x, phi_i);       % Integrale su tutta la lunghezza
        phi_at_xj(i) = phi_i(pos_xj);             % Valore nel punto osservato
        phiL(i) = phi_i(end);                     % Valore in x = L (estremo libero)
    end

    % Calcolo della risposta in frequenza teorica
    for k = 1:length(omega)
        w = omega(k);      % Frequenza attuale
        sum_modes_sy = 0;  % Inizializzazione somma modale

        for i = 1:length(i_nat)
            w_n = omega(i_nat(i));     % Frequenza naturale del modo i
            mi = modal_mass(i);        % Massa modale
            zeta_i = xsi;              % Smorzamento
            
            numerator = w^2 * d3_phi_at_0(i);  % Numeratore: w^2 * φ'''(0)
            denominator = mi * (-w^2 + 2j*zeta_i*w*w_n + w_n^2);  % Denominatore

            % Fattore di partecipazione
            participation = density * A * phi_integrals(i) + M_a * phiL(i);
            
            % Accumulo
            sum_modes_sy = sum_modes_sy + (numerator / denominator * participation);
        end

        G_sy(k) = -2 * E * J * sum_modes_sy;  % Somma finale
    end

    % Salva nella struttura G
    G.Gsy = G_sy;                   % Funzione complessa
    G.Gsy_amp = abs(G_sy);         % Modulo
    G.Gsy_phase = angle(G_sy);     % Fase
end

