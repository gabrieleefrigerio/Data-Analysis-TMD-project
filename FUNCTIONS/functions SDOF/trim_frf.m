clear
clc

% Finestra di dialogo per scegliere il file .mat
[filename, pathname] = uigetfile('data\*.mat', 'Seleziona un file MAT');
if isequal(filename, 0)
    disp('Nessun file selezionato.');
    return;
end

% Carica il file selezionato
load(fullfile(pathname, filename));

% Parametri di taglio
fcut = 480;
index = find(fvec < fcut);
frf = frf(index, :);
coh = coh(index, :);
fvec = fvec(index, :);

% Crea il nome del file di salvataggio aggiungendo "_trim"
[~, name, ~] = fileparts(filename); % Estrae il nome senza estensione
outputFilename = fullfile('Results', [name '_trim.mat']);

% Salva i risultati
save(outputFilename, "frf", "coh", "fvec", "fcut");
