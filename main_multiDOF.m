%% === MULTI-REFERENCE IBRAHIM TIME DOMAIN METHOD ===
% Author: Gabriele Frigerio

clc;
clear;
close all;

%% === STEP 1: Load FRF data ===
% Load Frequency Response Function (FRF) data from file
% Output:
%   - FRF: structure containing frequency-domain data
%   - acqNumbers: acquisition numbers
%   - filenames: names of the data files
[FRF, acqNumbers, filenames] = loadFRF_MDOF();

%% === STEP 2: Compute Impulse Response Function (IRF) ===
% Estimate the sampling frequency as twice the maximum frequency in FRF
fs = 2 * max(FRF.freq);

% Compute time-domain Impulse Response Function from FRF
[IRF] = computeIRF(FRF, fs);

%% === STEP 3: Detect the end of transient response ===
% Define threshold and window size for end-of-transient detection
fracThreshold = 0.0005;  % Relative threshold (e.g., 0.01% of peak)
windowSize    = 20;      % Number of consecutive samples for steady-state

% Automatically detect the index where the transient ends for each signal
kEnd = detectTransientEnd_minPerSignal(IRF, fracThreshold, windowSize);

%% === STEP 4: Build Stabilization Diagram and Perform Modal Analysis ===
% Extract IRF matrix dimensions: time samples x channels x acquisitions
h = IRF.irf;
[Nt, nChan, nAcq] = size(h);

% Maximum model order to consider
max_order = 20;

% Build the stabilization diagram and estimate modal parameters
modal_results = stabilization_diagram(IRF, kEnd, max_order, FRF);

%%

modal_results.modal_mass = find_modal_mass(FRF, modal_results);
freq = FRF.freq;
omega_vec = freq*2*pi;
FRF_fitted = zeros(length(omega_vec), 1)';
for ii = 1:length(modal_results.eigenfreq)
    FRF_fitted = FRF_fitted + (modal_results.modes(ii)*modal_results.modes(ii))./(modal_results.modal_mass(ii).*(-omega_vec.^2 + 2j*modal_results.damping(ii)*(modal_results.eigenfreq(ii)*2*pi).*omega_vec + (modal_results.eigenfreq(ii)*2*pi)^2));
end
figure;
semilogy(freq, abs(FRF_fitted), 'r--', freq, abs(FRF.FRF), 'b-');            
xlabel('Frequenza [Hz]');
ylabel('|FRF|');
legend('fitted FRF', 'experimental FRF');
grid on;


if ~exist('ibrahim_results', 'dir'), mkdir('ibrahim_results'); end; save('ibrahim_results/modal_results_I0.mat', 'modal_results');
