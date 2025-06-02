%% FITTING OF MODAL PARAMETERS

% Author: Gabriele Frigerio
% Last modified: 08/05/2025
%% Initialization
clc
clear
close all
set(0,'defaultaxesfontsize',16);
set(0, 'DefaultLineLineWidth', 1.5);

%% Load FRF data
[data, acqNumber, accIndex] = loadFRFData();
if isempty(data)
    return; % utente ha annullato la selezione
end

% Total number of channels
num_acc = size(data.frf, 2);

% Channel selection
accIndex = selectAccelerometer(num_acc, data);
%% Find FRF peaks
% [peaks, locs, locs_pos] = PeakSelectionGUI(f_range, FRF_range); % Avvia la GUI per la selezione dei picchi
peaks = [18 1.5]';
locs = [18 258];
locs_pos = [200 5000]';
%% Modes fitting

fitSingleMode(1 , peaks, locs, f_range, FRF_range, [], acqNumber, accIndex);

