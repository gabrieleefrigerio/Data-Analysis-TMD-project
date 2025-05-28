function [IRF] = computeIRF(data, fs)
    % computeIRF Trasforma la matrice FRF in IRF multi-acquisizione
    %   [IRF, time, freq] = computeIRF(data, fs)
    %   data: struct con campi
    %       - FRF: nFreq x nChan x nAcq
    %       - freq: nFreq x 1
    %   fs: frequenza di campionamento (Hz), opzionale

    % Verifica campi
    if ~isfield(data, 'FRF') || ~isfield(data, 'freq')
        error('Struct "data" deve contenere campi "FRF" e "freq".');
    end

    FRF = data.FRF;   % nFreq x nChan x nAcq
    freq = data.freq; % nFreq x 1

    % Stima fs se non fornito
    if nargin < 2 || isempty(fs)
        fs = 2 * max(freq);
    end

    [nFreq, nChan, nAcq] = size(FRF);
    Nfft = 2 * (nFreq - 1);  % lunghezza del segnale IRF

    % Preallocazione spettro completo
    FRF_full = zeros(Nfft, nChan, nAcq);
    % Parte positiva (0 → Nyquist)
    FRF_full(1:nFreq,:,:) = FRF;
    % Parte negativa (Nyquist→Fs)
    FRF_full(nFreq+1:end,:,:) = conj(flip(FRF(2:end-1,:,:), 1));

    % IFFT lungo la dimensione frequenza
    IRF = struct;
    IRF.irf = real(ifft(FRF_full, [], 1));   % Nfft x nChan x nAcq

    % Vettore tempo
    IRF.time = (0:Nfft-1)' / fs;
    IRF.freq = freq;
end