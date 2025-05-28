function [Phi, Qall] = scaling_modes(FRF, modal_results)
% scaling_modes - Scala i modi modali a partire da FRF e risultati modali
%
% Input:
%   FRF.freq       : [1 x M] frequenze (Hz)
%   FRF.FRF        : [M x nOut x nIn] funzione di risposta in frequenza
%   modal_results.modes : [nOut x 2N] modi non scalati
%   modal_results.poles : [1 x 2N] poli complessi
%
% Output:
%   Phi  : [nOut x 2N] modi scalati
%   Qall : {nOut x nIn} parametri stimati [Q; RU; RL] per ogni FRF

    % Estrai dati
    omega = 2 * pi * FRF.freq(:);  % [M x 1]
    H = FRF.FRF;                   % [M x nOut x nIn]
    psi = modal_results.modes;    % [nOut x 2N]
    s = modal_results.poles(:).'; % [1 x 2N]
    
    [M, nOut, nIn] = size(H);
    N2 = length(s);

    % Verifica dimensioni compatibili
    if size(psi, 2) ~= N2
        error('Il numero di colonne in psi deve corrispondere alla lunghezza di s.');
    end

    Qsum = zeros(1, N2);  % accumulatore per Q medi
    count = 0;
    Qall = cell(nOut, nIn);  % per salvare i Q stimati

    for p = 1:nOut
        for q = 1:nIn
            psi_p = psi(p, :);  % 1 x 2N
            psi_q = psi(q, :);  % <-- CORRETTO: ora prende q, non p

            H_pq = (H(:, p, q));  % [M x 1]
            A = zeros(M, N2 + 2);        % sistema lineare

            for i = 1:M
                jw = 1j * omega(i);
                for r = 1:N2
                    A(i, r) = (psi_p(r) * psi_q(r)) / (jw - s(r));
                end
                A(i, end-1) = 1;              % R_U
                A(i, end)   = 1 / omega(i)^2; % R_L
            end

            x = A \ H_pq;  % soluzione LS

            Qr = x(1:N2).';      % 1 x 2N
            Qall{p, q} = x;      % salva Q, R_U, R_L per questa FRF

            Qsum = Qsum + Qr;
            count = count + 1;
        end
    end

    Qavg = Qsum / count;         % media su tutte le FRF
    sqrtQ = sqrt(Qavg);          % 1 x 2N
    Phi = psi .* repmat(sqrtQ, size(psi, 1), 1);  % broadcasting esplicito
end
