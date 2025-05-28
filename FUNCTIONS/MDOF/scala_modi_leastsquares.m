function Psi_scaled = scala_modi_leastsquares(omega, Hpq_matrix, omega_r, xi, eta, Psi_unscaled)
% SCALA_MODI_LEASTSQUARES
%   Esegue la normalizzazione dei modi propri non scalati ottenuti da MITD,
%   basandosi su un fitting least-squares della FRF secondo la formulazione modale.

% Input:
%   omega         - vettore (k x 1) di frequenze [rad/s]
%   Hpq_matrix    - matrice (k x npairs) delle FRF misurate H_pq(omega) (una colonna per ogni accoppiamento p–q)
%   omega_r       - vettore (1 x N) delle frequenze naturali stimate
%   xi            - vettore (1 x N) degli smorzamenti stimati (smorzamento critico)
%   eta           - vettore (1 x N) dei fattori di scala temporali (dal MITD)
%   Psi_unscaled  - matrice (nDOF x N) dei modi propri non scalati (Ψ̂ = Ψ * Ψ̃)

% Output:
%   Psi_scaled    - matrice (nDOF x N) dei modi propri scalati

% --- Preliminari
[k, npairs] = size(Hpq_matrix);
N = length(omega_r);           % numero di modi
Psi_products = zeros(N, npairs);

% Loop su ogni accoppiamento p–q (una colonna di Hpq_matrix)
for pair = 1:npairs
    Hpq = Hpq_matrix(:, pair);   % vettore FRF per p–q

    % Costruzione dei termini modali nel dominio della frequenza
    D = zeros(k, N);
    for i = 1:N
        D(:, i) = 1 ./ (2j * eta(i) * omega_r(i) * sqrt(1 - xi(i)^2));
    end

    % Aggiungiamo componenti residui (costante e razionale)
    A = [D, ones(k,1), 1 ./ (omega.^2)];

    % Least squares fit per ottenere prodotti psi_p * psi_q
    x = A \ Hpq;

    % Salviamo i prodotti modali
    Psi_products(:, pair) = x(1:N);
end

% --- Ricostruzione dei modi propri scalati

% Supponiamo che ogni colonna in Psi_unscaled corrisponda ad un DoF
% e che Psi_products contenga i prodotti psi_p * psi_q
% Costruiamo una matrice dei prodotti Ψ_i Ψ_j^T
% Poi la decomponiamo per trovare Ψ scalati (fattorizzazione simmetrica)

% Stima numero di DoF da Psi_unscaled
nDOF = size(Psi_unscaled, 1);

% Ricostruiamo la matrice dei residui simmetrica (N x N)
R = zeros(N, N);

% Supponiamo che gli accoppiamenti siano p = q = 1, 2, ..., nDOF
% e ordinati come (1,1), (2,2), ..., (nDOF,nDOF)
% In alternativa, dovresti avere la lista degli indici (p,q) usati
count = 1;
for i = 1:nDOF
    for j = i:nDOF
        R(:, :) = R(:, :) + Psi_products(:, count) * (Psi_unscaled(i,:)' .* Psi_unscaled(j,:)');
        count = count + 1;
    end
end

% Decomposizione simmetrica (fattorizzazione Cholesky/SVD)
% per ottenere Ψ scalati
[U, S, ~] = svd(R);
Psi_scaled = Psi_unscaled * U * sqrt(S);

end
