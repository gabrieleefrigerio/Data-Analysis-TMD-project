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
fracThreshold = 0.0001;  % Relative threshold (e.g., 0.01% of peak)
windowSize    = 20;      % Number of consecutive samples for steady-state

% Automatically detect the index where the transient ends for each signal
kEnd = detectTransientEnd_minPerSignal(IRF, fracThreshold, windowSize);

%% === STEP 4: Build Stabilization Diagram and Perform Modal Analysis ===
% Extract IRF matrix dimensions: time samples x channels x acquisitions
h = IRF.irf;
[Nt, nChan, nAcq] = size(h);

% Maximum model order to consider
max_order = 25;

% Build the stabilization diagram and estimate modal parameters
modal_results = stabilization_diagram(IRF, kEnd, max_order, FRF);

%% === STEP 5: Normalize Mode Shapes ===

%% === or Scale Mode Shapes Using modal mass ===
% Perform scaling of mode shapes using FRF data
[Phi, Qall] = scaling_modes(FRF, modal_results);%% DA SISTEMARE
