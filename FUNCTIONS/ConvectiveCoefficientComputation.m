function [data] = ConvectiveCoefficientComputation(data, G)
%CONVECTIVECOEFFICIENTCOMPUTATION Calcola il coefficiente di scambio termico h
% basato su un flusso esterno forzato su un cilindro vibrante in aria.
%
% INPUT:
%   - data: struct contenente le proprietà del materiale e dell'aria
%   - G: funzione di trasferimento (FRF) complessa W(jΩ)/Y(jΩ)
%
% OUTPUT:
%   - data.h: coefficiente di scambio convettivo [W/(m^2·K)]

    % Diametro idraulico
    D_h = data.D_ext - data.D_int;  % [m]

    % === Proprietà dell'aria ===
    rho_a  = data.air.rho_a;     % [kg/m^3]
    mu     = data.air.ba;        % [Pa·s]
    cp     = data.air.eta_a;     % [J/(kg·K)]
    k_air  = data.air.fi_a;      % [W/(m·K)]

    % === Stima della velocità caratteristica ===
    % Assume G rappresenti FRF = W(jΩ)/Y(jΩ)
    % v = d(w)/dt → in freq: v = jω * G * Y
    % Qui stimiamo vc = media(|jω * G|) su tutte le frequenze

    omega = G.omega;        % rad/s (ipotesi frequenze analizzate)
    
    load("White_Noise_4_450_Hz.mat");

    
    Y_fft = fft(Signal_filtred);           % FFT di y(t)
    

    % Calcolo v in frequenza
    V_fft = 1i * omega .* G.Gwy .* Y_fft;

    % Passaggio nel dominio del tempo
    v_xt_T = ifft(V_fft);    % velocità v(x,t,T)
    % v_xt_T = real(ifft(V_fft)); oppure uso questo

    % Definizione dominio spazio-temporale
    L = data.L;               % lunghezza [m]
    
    load("time.mat");
    t_max = t(end);           % durata [s]
    Nx = 10000;              % punti spaziali
    Nt = length(t);             % punti temporali
    
    x = linspace(0, L, Nx);
    t = linspace(0, t_max, Nt);
    dx = x(2) - x(1);
    dt = t(2) - t(1);

    % Media temporale in ogni punto spaziale
    v_x_avgT = trapz(v_xt_T, 2) * dt / Tmax;    % risultato: Nx x 1

    % Media spaziale → v(T)
    v_T = trapz(v_x_avgT) * dx / L;
    % scalare
    v_c = v_T;

    % === Numeri adimensionali ===
    Re = rho_a * v_c * D_h / mu;                  % Reynolds number
    Pr = cp * mu / k_air;                       % Prandtl number

    % === Nusselt number con la correlazione per cilindro trasversale ===
    Nu = 0.3 + ...
        (0.62 * Re^(0.5) * Pr^(1/3)) / ...
        (1 + (0.4 / Pr)^(2/3))^(1/4) * ...
        (1 + (Re / 282000)^(5/8))^(4/5);

    % === Calcolo del coefficiente di scambio convettivo h ===
    h = Nu * k_air / D_h;                         % [W/(m^2·K)]

    % === Salvataggio ===
    data.h = h;
end



