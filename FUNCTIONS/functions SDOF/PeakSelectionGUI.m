function [peaks, locs, locs_pos] = PeakSelectionGUI(f_range, FRF_range)
%RUNPEAKSELECTIONGUI GUI per la selezione dei picchi su una FRF
%   Input:
%       - f_range: vettore delle frequenze
%       - FRF_range: vettore FRF (complesso o reale)
%   Output:
%       - peaks: valori dei picchi
%       - locs: posizioni in frequenza dei picchi
%       - locs_pos: indici nei dati originali

    peaks = [];
    locs = [];
    locs_pos = [];

    fig = uifigure('Name', 'Peak Finder', 'Position', [100 100 900 550]);

    uilabel(fig, 'Text', 'Ampiezza minima del picco:', 'Position', [30, 430, 200, 22]);
    defaultMinAmp = round(0.15 * max(abs(FRF_range)));
    minAmpField = uieditfield(fig, 'numeric', 'Value', defaultMinAmp, 'Position', [30, 400, 150, 22]);

    btnFind = uibutton(fig, 'Text', 'Trova Picchi', 'Position', [30, 360, 150, 30]);
    btnProceedToFit = uibutton(fig, 'Text', 'Procedi', 'Position', [30, 310, 150, 30]);

    ax = uiaxes(fig, 'Position', [200, 50, 670, 450]);
    title(ax, 'FRF con Picchi');
    xlabel(ax, 'Frequenza (Hz)');
    ylabel(ax, '|FRF|');
    grid(ax, 'on'); grid(ax, 'minor');

    btnFind.ButtonPushedFcn = @(~,~) findPeaksCallback();
    btnProceedToFit.ButtonPushedFcn = @(~,~) proceedCallback();

    % Lancio iniziale
    findPeaksCallback();
    uiwait(fig);

    %% === Callback Trova Picchi ===
    function findPeaksCallback()
        minAmp = minAmpField.Value;
        cla(ax);
        semilogy(ax, f_range, abs(FRF_range)); hold(ax, 'on');
        [pks, locsFound] = findpeaks(abs(FRF_range), f_range, 'MinPeakHeight', minAmp);
        [~, locsPos] = findpeaks(abs(FRF_range), 'MinPeakHeight', minAmp);
        semilogy(ax, locsFound, pks, 'ro', 'MarkerFaceColor', 'r');
        peaks = pks;
        locs = locsFound;
        locs_pos = locsPos;
        legend(ax, 'FRF', 'Picchi trovati');
        hold(ax, 'off');
    end

    function proceedCallback()
        if isempty(peaks)
            uialert(fig, 'Trova prima i picchi!', 'Errore');
            return;
        end
    
           uiresume(fig);  % sblocca il main
            close(fig);     % chiudi la GUI
    end

end
