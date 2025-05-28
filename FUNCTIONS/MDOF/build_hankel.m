function [H0, H1] = build_hankel(IRF, m, n)
% build_hankel - Costruisce matrici Hankel H0 e H1 secondo MRITD
%   [H0, H1] = build_hankel(IRF, m, n)
%
%   IRF: struttura con campo irf (Nt x nChan x nAcq)
%   m: numero di righe di blocchi temporali
%   n: numero di colonne di blocchi temporali
%   H0: Hankel a tempo t
%   H1: Hankel a tempo t + Δt (shiftata di 1)

    h = IRF.irf;  % Nt x nChan x nAcq
    [Nt, nChan, nAcq] = size(h);
    
    % Verifica sufficienza dati temporali
    maxShift = m + n - 1;
    if maxShift + 1 > Nt
        error('IRF troppo corta: richiede almeno %d campioni temporali.', maxShift + 1);
    end

    % Nuove dimensioni desiderate: m*nChan righe, n*nAcq colonne
    row_dim = m * nChan;
    col_dim = n * nAcq;

    % Preallocazione
    H0 = zeros(row_dim, col_dim);
    H1 = zeros(row_dim, col_dim);

    % Costruzione per blocchi: ogni blocco h(t + i) è una matrice nChan x nAcq
    for row = 1:m
        for col = 1:n
            time_idx0 = row + col - 1;      % tempo per H0
            time_idx1 = row + col;          % tempo per H1 (shiftato)

            % h(t): slice [nChan x nAcq]
            h0_blk = squeeze(h(time_idx0, :, :));  % [nChan x nAcq]
            h1_blk = squeeze(h(time_idx1, :, :));  % [nChan x nAcq]

            % Calcola indici di riga e colonna nel risultato
            r_start = (row - 1) * nChan + 1;
            r_end   = row * nChan;
            c_start = (col - 1) * nAcq + 1;
            c_end   = col * nAcq;

            % Inserisci i blocchi
            H0(r_start:r_end, c_start:c_end) = h0_blk;
            H1(r_start:r_end, c_start:c_end) = h1_blk;
        end
    end
end
