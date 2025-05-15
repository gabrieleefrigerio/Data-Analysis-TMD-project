function [phi_at_all_x, G] = GwyMultiComputation(omega, i_nat, modes_shapes, x, modal_mass, xsi, density, A, M_a, G)
% GWYCOMPUTATION Calcola la funzione di trasferimento teorica G_WY (acceleration FRF)
% per tutti i punti spaziali del vettore x.
%
% INPUT:
%   omega        - Vettore delle frequenze angolari (rad/s)
%   i_nat        - Indici delle frequenze naturali
%   modes_shapes - Matrice delle forme modali (una riga per ogni modo)
%   x            - Coordinate spaziali della trave
%   modal_mass   - Vettore delle masse modali
%   xsi          - Smorzamento modale (assunto costante per tutti i modi)
%   density      - Densità del materiale
%   A            - Area della sezione trasversale
%   M_a          - Massa concentrata all’estremo libero
%   G            - Struttura di output per salvare i risultati
%
% OUTPUT:
%   phi_at_all_x - Matrice (n_modi x n_punti_x) con le forme modali in ogni punto x
%   G            - Struttura con:
%                  Gwy: matrice (n_frequenze x n_punti_x) FRF complesse
%                  Gwy_amp: ampiezze (modulo)
%                  Gwy_phase: fasi (argomento)

    n_modes = length(i_nat);
    n_points = length(x);
    n_freq = length(omega);

    % Preallocazione
    G_FRF_theory = zeros(n_freq, n_points);  % Matrice FRF: righe -> frequenze, colonne -> punti spaziali
    phi_integrals = zeros(1, n_modes);       % Integrali spaziali delle forme modali
    phi_at_all_x = modes_shapes;            % Ogni riga è un modo, ogni colonna è un punto
    phiL = zeros(1, n_modes);                % Valore della forma modale all’estremo libero (x = L)

    % Calcolo degli integrali e del valore in x = L
    for i = 1:n_modes
        phi_i = modes_shapes(i, :);
        phi_integrals(i) = trapz(x, phi_i);   % Integrale forma modale
        phiL(i) = phi_i(end);                 % Valore all’estremo libero
    end

    % Ciclo su tutte le frequenze
    for k = 1:n_freq
        w = omega(k);  % Frequenza attuale

        % Ciclo su ogni punto spaziale xj
        for j = 1:n_points
            sum_modes = 0;

            for i = 1:n_modes
                w_n = omega(i_nat(i));
                mi = modal_mass(i);
                zeta_i = xsi;

                phi_ij = phi_at_all_x(i, j);  % Valore della forma modale nel punto j
                numerator = w^2 * phi_ij;
                denominator = mi * (-w^2 + 2j * zeta_i * w * w_n + w_n^2);
                participation = density * A * phi_integrals(i) + M_a * phiL(i);

                sum_modes = sum_modes + (numerator / denominator) * participation;
            end

            G_FRF_theory(k, j) = 1 + sum_modes;
        end
    end

    % Salvataggio nella struttura G
    G.Gwy = G_FRF_theory;               % FRF complessa (frequenze x posizioni)
    G.Gwy_amp = abs(G_FRF_theory);      % Ampiezze
    G.Gwy_phase = angle(G_FRF_theory);  % Fasi
end


