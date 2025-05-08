clc
clear variables
close all

%% SMA TMD PROJECT MAIN

%% INPUT CORRENT
I = 0:0.5:7.5; % input current for Transfer function computation

%% STRUCT INIZIALIZATION WITH PAPER PARAMETERS
data = ParametersInizialization();

%% RELATIONSHIP: Temperature and concentration
data = ConcentrationComputation(data);

%% COMPUTATION: Young Modulus and resistivity and damping
data.rho = data.ratio * data.rho_m + (1 - data.ratio) * data.rho_a;  % Resistività effettiva [Ohm·m]
data.E   = data.ratio * data.Em + (1 - data.ratio) * data.Ea;    % Modulo elastico effettivo [GPa]
data.xi1 =  data.ratio * data.xi1_m + (1 - data.ratio) * data.xi1_a;    % Coefficiente di damping adimensionale effettivo []

%% TRANSFER FUNCTION COMPUTATION
G = TransferFunctionComputation(data);

%% CALCULATION: h coefficient
data  = ConvectiveCoefficientComputation(data, G);

%% RELATIONSHIP: i Current and Temperature
data = PowerBalance(data);

%% CALCULATION: Transfer Function in dependence of i (current)

for ii = 1:length(I)
    % estabilished current value
    data.Current = Current(ii);

    % calculate temperature value
    data = PowerBalance(data);

    % calculate Martensite concentration
    data = ConcentrationComputation(data);

    % computation: Young Modulus and resistivity and damping
    data.rho = data.ratio * data.rho_m + (1 - data.ratio) * data.rho_a;  % Resistività effettiva [Ohm·m]
    data.E   = data.ratio * data.Em     + (1 - data.ratio) * data.Ea;    % Modulo elastico effettivo [GPa]
    data.xsi1 =  data.ratio * data.xi1_m + (1 - data.ratio) * data.xi1_a;    % Coefficiente di damping adimensionale effettivo []


    % TRANSFER FUNCTION COMPUTATION
    G = TransferFunctionComputation(data);

end

%% COMPUTATION: Shaker Signal

%% PLOT Transfer Function

%% PLOT Shaker Signal

%% EXPORT Transfer Function

%% EXPORT Shaker Signal 





