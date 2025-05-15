function [phi_at_xj, G] = GwySingleComputation(omega, i_nat, modes_shapes, x, pos_xj, modal_mass, xsi, density, A, M_a, G)
% GWYCOMPUTATION Calcola la funzione di trasferimento teorica G_WY (acceleration FRF).
%
% Questa funzione stima la risposta in accelerazione (funzione di trasferimento teorica) 
% in una posizione della trave (pos_xj), dovuta all'eccitazione modale. Utilizza il metodo 
% modale per approssimare la funzione di risposta in frequenza (FRF) in un punto.
%
% INPUT:
%   omega        - Vettore delle frequenze angolari (rad/s)
%   i_nat        - Indici delle frequenze naturali
%   modes_shapes - Matrice delle forme modali (una riga per ogni modo)
%   x            - Coordinate spaziali della trave
%   pos_xj       - Indice del punto di osservazione per la risposta (sull'asse x)
%   modal_mass   - Vettore delle masse modali
%   xsi          - Smorzamento modale (assunto costante per tutti i modi)
%   density      - Densità del materiale
%   A            - Area della sezione trasversale
%   M_a          - Massa concentrata all’estremo libero
%   G            - Struttura di output per salvare i risultati
%
% OUTPUT:
%   phi_at_xj    - Valori delle forme modali nel punto pos_xj
%   G            - Struttura con: Gwy (FRF complessa), Gwy_amp (ampiezza), Gwy_phase (fase)

    % Preallocazione dei risultati
    G_FRF_theory = zeros(length(omega), 1);        % FRF calcolata
    phi_integrals = zeros(1, length(i_nat));       % Integrali spaziali delle forme modali
    phi_at_xj = zeros(1, length(i_nat));           % Valore forma modale nel punto osservato
    phiL = zeros(1, length(i_nat));                % Valore forma modale all'estremo libero

    % Calcolo dei coefficienti modali necessari (integrali, valori locali)
    for i = 1:length(i_nat)
        phi_i = modes_shapes(i, :);
        phi_integrals(i) = trapz(x, phi_i);            % Integrale della forma modale
        phi_at_xj(i) = phi_i(pos_xj);                  % Valore nel punto osservato
        phiL(i) = phi_i(end);                          % Valore in x = L (estremo libero)
    end

    % Calcolo della risposta in frequenza teorica per ogni omega
    for k = 1:length(omega)
        w = omega(k);      % Frequenza attuale
        sum_modes = 0;     % Inizializzazione somma modale

        for i = 1:length(i_nat)
            w_n = omega(i_nat(i));      % Frequenza naturale del modo i
            mi = modal_mass(i);         % Massa modale del modo i
            zeta_i = xsi;               % Smorzamento modale (costante)

            % Numeratore: contributo del modo alla risposta in xj
            numerator = w^2 * phi_at_xj(i);

            % Denominatore: modello armonico smorzato per il modo i
            denominator = mi * (-w^2 + 2j * zeta_i * w * w_n + w_n^2);

            % Fattore di partecipazione del modo i
            participation = density * A * phi_integrals(i) + M_a * phiL(i);

            % Somma modale (superposizione lineare)
            sum_modes = sum_modes + (numerator / denominator) * participation;
        end

        % Funzione di risposta teorica
        G_FRF_theory(k) = 1 + sum_modes;
    end

    % Salvataggio nella struttura G
    G.Gwy = G_FRF_theory;                % FRF teorica (complessa)
    G.Gwy_amp = abs(G_FRF_theory);       % Ampiezza (modulo)
    G.Gwy_phase = angle(G_FRF_theory);   % Fase (argomento)
end