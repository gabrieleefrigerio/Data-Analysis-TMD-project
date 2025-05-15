function [norm_sp, freq_vec]=fft_n(data,fsamp)
% Questa funzione calcola la FFT normalizzata di un segnale (o più segnali)
% e restituisce solo le frequenze positive.
%
% INPUT:
%   - data  : matrice dei dati [r campioni, c segnali]
%   - fsamp : frequenza di campionamento [Hz]
%
% OUTPUT:
%   - norm_sp : spettro normalizzato (solo frequenze positive)
%   - freq_vec: vettore delle frequenze corrispondenti alle componenti spettrali

    % Determina la dimensione della matrice dei dati
    dim = size(data);

    % Verifica che i segnali siano organizzati in colonna.
    % Se sono in riga (più colonne che righe), li trasforma in colonna.
    if dim(2) > dim(1)
        data = data';
    end

    % Numero totale di campioni nel segnale
    N = length(data);

    % Risoluzione spettrale, dipendente dalla lunghezza del segnale
    df = fsamp / N;

    % Verifica se il numero di campioni è pari o dispari
    if (N/2) == floor(N/2)
        % --- Caso N pari ---

        % Crea il vettore delle frequenze da 0 fino a fs/2
        freq_vec = [0:df:(N/2 * df)]';

        % Numero di frequenze positive
        NF = length(freq_vec);

        % Calcola la FFT lungo la prima dimension (cioè per colonna)
        sp = fft(data, [], 1);

        % Normalizza il primo valore (DC) sul numero totale di campioni
        norm_sp(1,:) = sp(1,:) / N;

        % Normalizza i valori dal 2° fino al valore a fs/2 su metà del numero di campioni
        norm_sp(2:N/2,:) = sp(2:N/2,:) / (N/2);

        % Normalizza il valore a fs/2 (Nyquist) sul numero totale di campioni
        norm_sp(N/2+1,:) = sp(N/2+1,:) / N;

    else
        % --- Caso N dispari ---

        % Crea il vettore delle frequenze da 0 fino a (N-1)/2 * df
        freq_vec = [0:df:((N-1)/2) * df]';

        % Numero di frequenze positive
        NF = length(freq_vec);

        % Calcola la FFT lungo la prima dimension (cioè per colonna)
        sp = fft(data, [], 1);

        % Normalizza il primo valore (DC) sul numero totale di campioni
        norm_sp(1,:) = sp(1,:) / N;

        % Normalizza i valori dal 2° fino a (N+1)/2 su metà del numero di campioni
        norm_sp(2:(N+1)/2,:) = sp(2:(N+1)/2,:) / (N/2);
    end

end