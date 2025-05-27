function [data, acqNumber, accIndex] = loadFRFData()
    %% Load FRF data
    [filename, pathname] = uigetfile('*.mat', 'Seleziona un file MAT');
    if isequal(filename, 0)
        disp('Selezione file annullata.');
        data = [];
        acqNumber = [];
        accIndex = [];
        return;
    end

    % Carica il file
    data = load(fullfile(pathname, filename));

    % Estrai il numero acquisizione dal nome del file
    acqMatch = regexp(filename, 'acquisition_(\d+)', 'tokens');
    if ~isempty(acqMatch)
        acqNumber = str2double(acqMatch{1}{1});
    else
        warning('Numero acquisizione non trovato nel nome del file. Imposto a 1.');
        acqNumber = 1;  % fallback
    end

    % Estrai numero massimo accelerometri
    num_acc = size(data.frf, 2);

    % Selezione accelerometro
    accIndex = selectAccelerometer(num_acc, data);
end
