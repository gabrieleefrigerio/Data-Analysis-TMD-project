function [FRF, filename, pathname] = loadFRF_MDOF()
    % Load a single .mat file containing FRF data
    % Returns:
    %   - FRF: structure with FRF.H1, FRF.freq, FRF.gamma2
    %   - filename: name of the loaded file
    %   - pathname: folder containing the file

    [filename, pathname] = uigetfile('*.mat', 'Select a MAT file with FRF');
    
    if isequal(filename, 0)
        disp('File selection canceled.');
        FRF = [];
        return;
    end
    
    fileData = load(fullfile(pathname, filename));
    
    % Check for required variables
    if ~isfield(fileData, 'H1') || ~isfield(fileData, 'frq')
        error('The selected file does not contain required variables: "H1" and "frq".');
    end

    % Optional check for coherence
    if isfield(fileData, 'gamma2')
        gamma2 = fileData.gamma2;
    else
        gamma2 = [];
        warning('No coherence "gamma2" found in file.');
    end

    % Build output structure
    FRF.FRF     = fileData.H1;
    FRF.freq   = fileData.frq;
    FRF.gamma2 = gamma2;
end
