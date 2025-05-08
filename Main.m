clc
clear variables
close all

%% SMA TMD PROJECT MAIN

%% INPUT CORRENT
I = 0:0.5:7.5; % input current for Transfer function computation

%% PAPER PARAMETERS
% Creo la struct con tutti i dati necessari (presi dal paper) per le
% simulazioni
data = ParametersInizialization();
% Corrente all'interno della trave
data.Current = I(1);           % [A] corrente per cambiare la temperatura

%% RELATIONSHIP: Temperature and concentration
data = ConcentrationComputation(data);

%% COMPUTATION: Young Modulus and resistivity
data.rho = data.ratio * data.rho_m + (1 - data.ratio) * data.rho_a;  % Resistività effettiva [Ohm·m]
data.E   = data.ratio * data.Em     + (1 - data.ratio) * data.Ea;    % Modulo elastico effettivo [GPa]

%% TRANSFER FUNCTION COMPUTATION
[G_free, G_forced] = TransferFunctionComputation(data);

%% CALCULATION: h coefficient
data  = ConvectiveCoefficientComputation(data, G);

%% RELATIONSHIP: i Current and Temperature
data = PowerBalance(data);

%% CALCULATION: Transfer Function in dependence of i (current)

for ii = 1:length(I)
    
    data.Current = Current(ii);

end

%% COMPUTATION: Shaker Signal

%% PLOT Transfer Function

%% PLOT Shaker Signal

%% EXPORT Transfer Function

%% EXPORT Shaker Signal 





