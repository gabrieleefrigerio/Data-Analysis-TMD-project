%% SMA TMD PROJECT MAIN
clc
clear variables
close all

%% INPUT CORRENT
% Current = 0:0.1:8; % input current for Transfer function computation
Current = [0 2 5 5.8:0.1:6.3 7];

%% STRUCT INIZIALIZATION WITH PAPER PARAMETERS
data = ParametersInizialization();

% [°C] Temperatura Tmd iniziale
data.T = 25;                


%% CALCULATION: Transfer Function in dependence of i (current)
flag = 0;
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
    
    % SALVO GLI OUTPUT IN UNA STRUCT
    Output.Current(ii) = Current(ii);
    Output.Temperature(ii) = data.T;
    Output.ratio(ii) = data.ratio;
    Output.NaturalFrequencies(ii,:) = G.freq_nat;


    % TROVO LO SHIFT DI CORRENTE
    if data.ratio < 1 && flag == 0
        Output.TransStart = Output.Current(ii);
        flag = 1;
    elseif data.ratio > 0 && flag == 1
        Output.TransFinish = Output.Current(ii);
    end

end


%% PLOT FRF
GPlot(G.freq, Current, GSY_values, GWY_values)

%% DISP OUTPUT

disp("------ Results ------");
disp(Output);

disp("------ Natural Frequency Start [Hz] ------ ");
disp(Output.NaturalFrequencies(1,:));
disp("------ Natural Frequency Finish [Hz] ------ ");
disp(Output.NaturalFrequencies(end,:));