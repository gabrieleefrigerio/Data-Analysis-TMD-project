%% SMA TMD PROJECT MAIN
clc
clear variables
close all



%% RATIOS VECTOR
ratios = [0 1];

%% STRUCT INIZIALIZATION WITH PAPER PARAMETERS
data = ParametersInizialization();

for ii =1:length(ratios)
    %% DEFINE: Martensite ratio
    data.ratio = ratios(ii);
    %% COMPUTATION: Young Modulus and resistivity and damping
    data.rho = data.ratio * data.rho_m + (1 - data.ratio) * data.rho_a;  % Resistività effettiva [Ohm·m]
    data.E   = data.ratio * data.Em + (1 - data.ratio) * data.Ea;    % Modulo elastico effettivo [GPa]
    data.xi1 =  data.ratio * data.xi1_m + (1 - data.ratio) * data.xi1_a;    % Coefficiente di damping adimensionale effettivo []
    
    %% TRANSFER FUNCTION COMPUTATION
    G = TransferFunctionComputation(data, true);
    %% CALCULATION: h coefficient
    data  = ConvectiveCoefficientComputation(data, G);

    if ratios(ii) == 1
        data.h_m = data.h; % coef scambio convettivo martensite
    else
        data.h_a = data.h; % coef scambio convettivo austenite
    end
end

%% DISP RISULTATI
disp(" ------- SCAMBIO CONVETTIVO MARTENSITE [W/(m^2*K)]-------");
disp(data.h_m);
disp(" ------- SCAMBIO CONVETTIVO AUSTENITE [W/(m^2*K)]-------");
disp(data.h_a);