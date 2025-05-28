function [autocross_mean, frequencies] = autocross(data1, data2, fsamp, N_win, N_OL, Win)
% AUTOCROSS - Calcola lo spettro incrociato (o autospettro) tra due segnali
%             usando segmenti sovrapposti e una finestra specificata
%
% INPUT:
%   data1   - Primo segnale (vettore colonna)
%   data2   - Secondo segnale (vettore colonna)
%   fsamp   - Frequenza di campionamento [Hz]
%   N_win   - Lunghezza della finestra (in campioni)
%   N_OL    - Numero di campioni di overlap tra finestre
%   Win     - Finestra da applicare ai dati (es. hanning(N_win))
%
% OUTPUT:
%   autocross_mean - Spettro medio (auto o incrociato)
%   frequencies    - Vettore delle frequenze corrispondenti [Hz]

% Numero di campioni totali nel segnale
N = length(data1);  % data1 e data2 devono avere stessa lunghezza

% Risoluzione in frequenza
df = fsamp / N_win;

% Costruzione del vettore delle frequenze (solo semispettro positivo)
if mod(N_win, 2) == 0
    frequencies = 0 : df : (N_win / 2) * df;
else
    frequencies = 0 : df : ((N_win - 1) / 2) * df;
end

NF = length(frequencies);  % Numero di frequenze considerate

% Numero di segmenti (finestre) da analizzare
num_records = fix((N - N_OL) / (N_win - N_OL));

% Preallocazione matrice spettri incrociati (una colonna per finestra)
autocross = zeros(NF, num_records);

% Inizializzazione del contatore
counter = 1;
finalPoint_nextIT = 0;

% Ciclo su ciascuna finestra (con overlap)
while finalPoint_nextIT <= N
    % Inizio e fine della finestra corrente
    start_p = (counter - 1) * (N_win - N_OL) + 1;
    finish_p = start_p + (N_win - 1);

    % Estrai segmento di dati e applica la finestra
    seg1 = Win .* data1(start_p : finish_p);
    seg2 = Win .* data2(start_p : finish_p);

    % Calcolo FFT e normalizzazione
    sp1 = fft(seg1) / N_win;
    sp2 = fft(seg2) / N_win;

    % Calcolo del prodotto spettrale (spettro incrociato)
    % Solo le frequenze positive
    autocross(:, counter) = conj(sp1(1:NF)) .* sp2(1:NF);

    % Aggiorna contatore e punto finale della prossima iterazione
    counter = counter + 1;
    finalPoint_nextIT = finish_p + N_win - N_OL;
end

% Correzione dell'ampiezza: raddoppia le componenti non-DC e non-Nyquist
if mod(N_win, 2) == 0
    autocross(2:end-1, :) = autocross(2:end-1, :) * 2;
else
    autocross(2:end, :) = autocross(2:end, :) * 2;
end

% Media tra tutte le finestre
autocross_mean = mean(autocross(1:NF, :), 2);

end
