function [data] = ConvectiveCoefficientComputation(data, G)
% CONVECTIVECOEFFICIENTCOMPUTATION Calcola il coefficiente di scambio termico convettivo (h)
% per un cilindro cavo soggetto a vibrazioni trasversali in aria.
%
% INPUT:
%   - data: struttura contenente le proprietà geometriche del cilindro e
%           le proprietà termofisiche dell'aria
%   - G: struttura contenente la funzione di risposta in frequenza (FRF)
%        G.Gwy = W(jω)/Y(jω), e il vettore delle frequenze angolari G.omega
%
% OUTPUT:
%   - data.h: coefficiente di scambio termico convettivo medio [W/(m^2·K)]

    %% === 1. Calcolo del diametro idraulico ===
    D_h = data.D_ext - data.D_int;  % Diametro idraulico [m]

    %% === 2. Estrazione delle proprietà dell'aria dal campo data.air ===
    rho_a  = data.air.rho_a;   % Densità dell'aria [kg/m^3]
    mu     = data.air.ba;      % Viscosità dinamica dell'aria [Pa·s]
    cp     = data.air.eta_a;   % Calore specifico a pressione costante [J/(kg·K)]
    k_air  = data.air.fi_a;    % Conducibilità termica [W/(m·K)]

    %% === 3. Stima della velocità caratteristica media v_c ===

    % Caricamento dei segnali di ingresso
    load("White_Noise_4_450_Hz.mat");  % Segnale temporale di spostamento filtrato (Signal_filtred)
    load("time.mat");                  % Vettore temporale t

    f_samp = 1 / (t(2) - t(1));        % Frequenza di campionamento [Hz]
    t_max = t(end);                    % Durata totale del segnale [s]

    % Calcolo FFT dello spostamento nel tempo Y(t) → Y(ω)
    [Y_fft, freq_vec] = fft_n(Signal_filtred, f_samp);  % FFT di y(t)


    % Calcolo della velocità nel dominio della frequenza:
    % v(ω) = jω * G(ω) * Y(ω)
    % V_fft = 1i * G.omega.' .* G.Gwy .* Y_fft;  % parte reale della velocità
    V_fft = 1i * (G.omega .* Y_fft) .* G.Gwy;
    
    % === Calcolo N ===
    N = 2 * (length(V_fft) - 1);
    
    % === Denormalizza ===
    V_half = V_fft;
    V_half(2:end-1,:) = V_half(2:end-1,:) /  2;
    
    % === Ricostruzione spettro completo ===
    V_full = [conj(flipud(V_half(2:end-1,:))); V_half];
    
    % === Frequenze ===
    fs = f_samp;  % metti il tuo valore reale se lo conosci
    f_fft = linspace(0, fs/2, length(V_fft));  % frequenze positive (unilaterale)
    f_full = linspace(-freq_vec(end), freq_vec(end), N);               % spettro completo
    
    % % === Plot ===
    % figure;
    % 
    % % Spettro completo
    % plot(f_full, abs(V_full), 'r-', 'DisplayName', '|V\_full| (completo)');
    % hold on;
    % 
    % % Spettro unilaterale (solo parte positiva)
    % plot(f_fft, abs(V_fft), 'b--o', 'DisplayName', '|V\_fft| (unilaterale, normalizzato)');
    % 
    % title('Modulo spettro: completo vs unilaterale');
    % xlabel('Frequenza (Hz)');
    % ylabel('Ampiezza');
    % legend;
    % grid on;

    %figure; semilogy(f_full, abs(V_full)); hold on; semilogy(freq_vec, abs(V_fft));

    % Trasformazione nel dominio del tempo → v(x,t,T)
    v_xt_T = real(ifft(V_full));  % Ricostruzione temporale della velocità [m/s]

    %% === 4. Definizione del dominio spazio-temporale ===
    L = data.L;                  % Lunghezza del cilindro [m]
    t_trapz = linspace(0,t(end), N);
    x = G.x;
    %% === 5. Calcolo della velocità media caratteristica v_c ===

    % Media temporale della velocità in ogni punto spaziale (integrale temporale)
    v_x_avgT = trapz(t_trapz,v_xt_T) / t_max;  % vettore Nx x 1

    % Media spaziale della velocità media → stima finale della velocità caratteristica
    v_T = trapz(x, v_x_avgT) / L;  % scalare
    v_c = v_T;

    %% === 6. Calcolo dei numeri adimensionali ===

    Re = rho_a * v_c * D_h / mu;    % Numero di Reynolds
    Pr = cp * mu / k_air;           % Numero di Prandtl

    %% === 7. Calcolo del numero di Nusselt (correlazione per cilindri trasversali) ===
    Nu = 0.3 + ...
         (0.62 * Re^0.5 * Pr^(1/3)) / ...
         (1 + (0.4 / Pr)^(2/3))^(1/4) * ...
         (1 + (Re / 282000)^(5/8))^(4/5);

    %% === 8. Calcolo del coefficiente convettivo h ===
    h = Nu * k_air / D_h;  % Coefficiente di scambio termico [W/(m^2·K)]

    %% === 9. Salvataggio del risultato ===
    data.h = h;  % Salva h nella struttura di output
end


