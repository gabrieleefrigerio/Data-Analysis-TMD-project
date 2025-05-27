function showFinalResult(allParams, f_range, FRF_range, accIndex, acqNumber)
    figFinal = figure('Name', 'FRF SDOF e Avvio MultiDOF', 'Position', [100 100 1600 900]);

    % Variabili locali
    omega_vec = 2 * pi * f_range;
    G_total = zeros(size(omega_vec));
    modes = cell(size(allParams, 1), 6);

    % Calcolo G_total e salvataggio parametri modali
    for i = 1:size(allParams, 1)
        p = allParams(i, :);
        A = p(3);
        G_total = G_total + A ./ (-omega_vec.^2 + 2j * p(2) * p(1) * omega_vec + p(1)^2);
        modes(i, :) = {p(1)/(2*pi), p(2), p(1), p(3), p(4), p(5)};
    end
    G_total = G_total + allParams(1, 5)./omega_vec.^2;

    % Ottimizza residui e aggiorna G_total
    [G_total, allParams, modes] = ottimizzaResidui(G_total, allParams, FRF_range, omega_vec, modes);

    % Tabella a sinistra
    modeTable = cell2table(modes, 'VariableNames', {'Frequenza (Hz)', 'Damping (xi)', 'omega0', 'A', 'Rh', 'Rl'});
    uitable(figFinal, 'Units', 'normalized', ...
        'Position', [0.05 0.15 0.20 0.75], ...
        'Data', modeTable{:,:}, ...
        'ColumnName', modeTable.Properties.VariableNames);

    % Grafico ampiezza
    axAmplitude = axes(figFinal, 'Units', 'normalized', 'Position', [0.35 0.55 0.6 0.4]);
    semilogy(axAmplitude, f_range, abs(FRF_range), 'b'); hold on;
    semilogy(axAmplitude, f_range, abs(G_total), 'r--', 'LineWidth', 1.5);
    xlabel(axAmplitude, 'Frequenza (Hz)'); ylabel(axAmplitude, 'Ampiezza');
    title(axAmplitude, 'Ampiezza');
    legend(axAmplitude, 'FRF', 'Somma SDOF'); grid(axAmplitude, 'on'); grid(axAmplitude, 'minor');

    % Grafico fase
    axPhase = axes(figFinal, 'Units', 'normalized', 'Position', [0.35 0.15 0.6 0.3]);
    plot(axPhase, f_range, angle(FRF_range), 'b'); hold on;
    plot(axPhase, f_range, angle(G_total), 'r--', 'LineWidth', 1.5);
    xlabel(axPhase, 'Frequenza (Hz)'); ylabel(axPhase, 'Fase (rad)');
    title(axPhase, 'Fase');
    legend(axPhase, 'FRF', 'Somma SDOF'); grid(axPhase, 'on'); grid(axPhase, 'minor');

    % Pulsanti
    uicontrol(figFinal, 'Style', 'pushbutton', 'String', 'Esporta Dati', ...
        'Position', [700 10 200 30], ...
        'Callback', @(~,~) exportData(G_total, f_range, modes, accIndex, acqNumber));

    uicontrol(figFinal, 'Style', 'pushbutton', 'String', 'Quit', ...
        'Position', [940 10 200 30], ...
        'Callback', @(~,~) close(figFinal));

    %% OTTIMIZZAZIONE RESIDUI FINALE
    function [G_total, allParams, modes] = ottimizzaResidui(G_total, allParams, FRF_range, omega_vec, modes)
        G_base = zeros(size(omega_vec));
        for i = 1:size(allParams,1)
            p = allParams(i,:);
            A = p(3);
            G_base = G_base + A ./ (-omega_vec.^2 + 2j*p(2)*p(1)*omega_vec + p(1)^2);
        end
        G_base = G_base + allParams(1,5)./omega_vec.^2;
        G_base = G_base.';
    
        r0 = allParams(end,4);
        opts = optimoptions('lsqnonlin', 'Display', 'iter', 'TolFun',1e-12, 'TolX',1e-12);
        obj = @(res) sum(real(G_base + res - FRF_range).^2 + imag(G_base + res - FRF_range).^2);
        res_ottimizzati = lsqnonlin(obj, r0, [], [], opts);
    
        allParams(end,4) = res_ottimizzati;
        modes(end,5) = {res_ottimizzati};
        G_total = G_base + res_ottimizzati;
    end

%% EXPORT DATA
    function exportData(G_total, f_range, modes, accIndex, acqNumber)
        % Percorso base
        resultsDirFRF = 'Results/Fitted FRFs';
        
        % Crea la cartella base se non esiste
        if ~exist(resultsDirFRF, 'dir')
            mkdir(resultsDirFRF);
        end
    
        % Crea la sottocartella per l'acquisizione
        acqFolderName = sprintf('Acquisizione_%d', acqNumber);
        resultsDirFRF = fullfile(resultsDirFRF, acqFolderName);
        
        if ~exist(resultsDirFRF, 'dir')
            mkdir(resultsDirFRF);
        end
    
        % Stessa cosa per la cartella delle tabelle
        resultsDirModes = fullfile(resultsDirFRF, 'Table Modi');
        if ~exist(resultsDirModes, 'dir')
            mkdir(resultsDirModes);
        end
    
        % Salva FRF
        suffix = sprintf('acc_%d', accIndex);
        frfMATname = fullfile(resultsDirFRF, ['FRF_SDOF_Optimize_' suffix '.mat']);
        frfData = G_total(:);
        freqData = f_range(:);
        save(frfMATname, 'frfData', 'freqData');
    
        % Salva Tabella Modi
        modeTable = cell2table(modes, ...
            'VariableNames', {'Frequenza_Hz', 'Damping_xi', 'omega0', 'A', 'Rh', 'Rl'});
        modesMATname = fullfile(resultsDirModes, ['ModeTable_FRF_' suffix '.mat']);
        save(modesMATname, 'modeTable');
        save(fullfile(resultsDirModes, 'freq'), 'f_range');
    
        msgbox({'Dati esportati con successo!'}, 'Esportazione completata');
    end

end
