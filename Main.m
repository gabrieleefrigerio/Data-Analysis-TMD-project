%% SMA TMD PROJECT MAIN
clc
clear variables
close all

%% INPUT CORRENT
Current = 0:0.5:7.5; % input current for Transfer function computation


%% STRUCT INIZIALIZATION WITH PAPER PARAMETERS
data = ParametersInizialization();

data.T = 25;                % [°C] Temperatura Tmd iniziale


%% CALCULATION: Transfer Function in dependence of i (current)

for ii = 1:length(Current)
    if ii == 1
        % calculate Martensite concentration
        data = ConcentrationComputation(data);
    
        % estabilished current value
        data.Current = Current(ii);
    
        % trovo il coef di scambio convettivo
        data.h = data.ratio * data.h_m + (1 - data.ratio) * data.h_a;  
        % computation: Young Modulus and resistivity and damping
        data.rho = data.ratio * data.rho_m + (1 - data.ratio) * data.rho_a;  % Resistività effettiva [Ohm·m]
        data.E   = data.ratio * data.Em     + (1 - data.ratio) * data.Ea;    % Modulo elastico effettivo [GPa]
        data.xi1 =  data.ratio * data.xi1_m + (1 - data.ratio) * data.xi1_a;    % Coefficiente di damping adimensionale effettivo []
        
        % calculate temperature value
        data = PowerBalance(data);
            
    else
        % estabilished current value
        data.Current = Current(ii);

        % trovo il coef di scambio convettivo
        data.h = data.ratio * data.h_m + (1 - data.ratio) * data.h_a;  

        % calculate temperature value
        data = PowerBalance(data);

        % calculate Martensite concentration
        data = ConcentrationComputation(data);

        % computation: Young Modulus and resistivity and damping
        data.rho = data.ratio * data.rho_m + (1 - data.ratio) * data.rho_a;  % Resistività effettiva [Ohm·m]
        data.E   = data.ratio * data.Em     + (1 - data.ratio) * data.Ea;    % Modulo elastico effettivo [GPa]
        data.xi1 =  data.ratio * data.xi1_m + (1 - data.ratio) * data.xi1_a;    % Coefficiente di damping adimensionale effettivo []
    end

    % TRANSFER FUNCTION COMPUTATION
    G = TransferFunctionComputation(data, false);

    % MEMORIZZA I VALORI DI G_SY E G_WY
    GSY_values(:,ii) = G.Gsy;   % Salva il valore di G_SY per la corrente corrente
    GWY_values(:,ii) = G.Gwy;   % Salva il valore di G_WY per la corrente corrente

    
end


%% PLOT FRF
GPlot(G.freq, Current, GSY_values, GWY_values)
