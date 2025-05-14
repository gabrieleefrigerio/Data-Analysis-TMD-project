function GPlotSingle(freq, G)
% GPLOTSINGLE Crea i grafici separati di ampiezza e fase per le funzioni di risposta G_WY e G_SY.
%
% Questa funzione disegna quattro grafici distinti:
%   1. Ampiezza della funzione teorica G_WY (accelerazione)
%   2. Fase della funzione teorica G_WY (in gradi)
%   3. Ampiezza della funzione sperimentale G_SY (accelerazione)
%   4. Fase della funzione sperimentale G_SY (in gradi)
%
% INPUT:
%   freq - Vettore delle frequenze [Hz]
%   G    - Struttura contenente:
%            G.Gwy     -> funzione teorica complessa (WY)
%            G.Gsy     -> funzione sperimentale complessa (SY)
%
% NOTA: Le funzioni devono essere già calcolate prima della chiamata.

    % === G_WY - Ampiezza ===
    figure('Color','w','WindowStyle','docked');  % Crea nuova figura ancorata
    semilogy(freq, abs(G.Gwy), 'b', 'LineWidth', 1.6);  % Traccia il modulo (ampiezza)
    grid on;
    xlabel('Frequenza [Hz]');
    ylabel('|G_{WY}|');
    title('Ampiezza G_{WY}');

    % === G_WY - Fase ===
    figure('Color','w','WindowStyle','docked');
    plot(freq, unwrap(angle(G.Gwy)) * 180/pi, 'b', 'LineWidth', 1.6);  % Fase in gradi (unwrap evita salti di 360°)
    grid on;
    xlabel('Frequenza [Hz]');
    ylabel('Fase G_{WY} [°]');
    title('Fase G_{WY}');

    % === G_SY - Ampiezza ===
    figure('Color','w','WindowStyle','docked');
    semilogy(freq, abs(G.Gsy), 'r', 'LineWidth', 1.6);  % Modulo sperimentale
    grid on;
    xlabel('Frequenza [Hz]');
    ylabel('|G_{SY}|');
    title('Ampiezza G_{SY}');

    % === G_SY - Fase ===
    figure('Color','w','WindowStyle','docked');
    plot(freq, unwrap(angle(G.Gsy)) * 180/pi, 'r', 'LineWidth', 1.6);  % Fase sperimentale
    grid on;
    xlabel('Frequenza [Hz]');
    ylabel('Fase G_{SY} [°]');
    title('Fase G_{SY}');
end
