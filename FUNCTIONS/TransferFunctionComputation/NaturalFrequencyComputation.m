function [i_nat, G] = NaturalFrequencyComputation(dets, freq)
% NATURALFREQUENCYCOMPUTATION Trova le frequenze naturali del sistema.
% Identifica gli indici corrispondenti ai minimi locali del determinante
% della matrice H (che indicano possibili modi propri).
%
% INPUT:
%   dets  - Vettore dei determinanti della matrice H(ω)
%   freq  - Vettore delle frequenze corrispondenti [Hz]
%
% OUTPUT:
%   i_nat - Indici delle frequenze naturali (minimi locali del det)
%   G     - Struttura contenente le frequenze naturali individuate

    % Inizializza vettore degli indici delle frequenze naturali
    i_nat = [];

    % Scansione dei valori del determinante per trovare i minimi locali
    for i = 2:length(dets)-1
        % Condizione per minimo locale (valore assoluto)
        if abs(dets(i)) < abs(dets(i-1)) && abs(dets(i)) < abs(dets(i+1))
            i_nat(end+1) = i;  % Aggiunge l'indice alla lista
        end
    end

    % Estrazione delle frequenze naturali dai relativi indici
    G.freq_nat = freq(i_nat);

    % Output delle frequenze naturali trovate
    % fprintf('Frequenze naturali [Hz]:\n');
    % disp(freq(i_nat));

    % --- (Facoltativo) Plot del determinante e delle frequenze naturali ---
    %{
    figure, box on;
    semilogy(freq, abs(dets), '-b');  % Plotta il determinante in scala logaritmica
    hold on;
    grid on;
    grid minor;
    ylabel('|det(H)|');
    xlabel('Frequenza [Hz]');
    title('Determinante della matrice H');
    plot(freq(i_nat), abs(dets(i_nat)), 'or');  % Evidenzia i minimi locali
    %}
end


