clc                 % Pulisce la finestra comandi
clear               % Cancella tutte le variabili dall'ambiente
close all           % Chiude tutte le figure aperte

%% ===================== IMPORT DATI =====================

data = struct();    % Inizializza una struttura vuota per contenere i dati caricati

% Caricamento dei vari file .mat in campi della struct 'data'
data.I0   = load('Prova8_6_4min\test@3830835382.mat');          
data.I1   = load('Corrente_1A\Prova1_4min\test@3830837208.mat');
data.I2   = load('Corrente_2A\Prova1_4min\test@3830837967.mat');
data.I3   = load('Corrente_3A\Prova1_4min\test@3830838831.mat');
data.I4   = load('Corrente_4A\Prova1_4min\test@3830839736.mat');
data.I5   = load('Corrente_5A\Prova1_4min\test@3830840454.mat');
data.I5_5 = load('Corrente_5A\Prova1_4min\test@3830840454.mat');
data.I5_8 = load('Corrente_5.8A\Prova1_4min\test@3830841870.mat');
data.I6_1 = load('Corrente_6.1A\Prova1_4min\test@3830842924.mat');
data.I6_4 = load('Corrente_6.4A\Prova2_4min\test@3830843986.mat');
data.I6_7 = load('Corrente_6.7A\Prova2_4min\test@3830844967.mat');
data.I7   = load('Corrente_7A\Prova2_4min\test@3830845894.mat');
data.I7_6 = load('Corrente_7.6A\Prova2_4min\test@3830847646.mat');
data.I8   = load('Corrente_8A\Prova2_4min\test@3830848913.mat');
data.I8_5 = load('Corrente_8.5A\Prova2_4min\test@3830849761.mat');
data.I9   = load('Corrente_9A\Prova1_4min\test@3830850372.mat');

%% ============= PARAMETRI ANALISI SPETTRALE =============

f_start_plot = 5;   % Frequenza minima per i plot
Time = 20;          % Durata della finestra di analisi in secondi
ol = 0.66;          % Sovrapposizione percentuale tra finestre

sens_force = 22.4e-3;   % Sensibilità sensore di forza [V/N]
sens_acc = 10.2e-3;     % Sensibilità accelerometro [V/(m/s^2)]

fsamp = 3125;           % Frequenza di campionamento [Hz]
dt = 1/fsamp;           % Intervallo temporale tra campioni [s]

%% ========== CALCOLO FRF E SALVATAGGIO IN STRUTTURA ==========

H1_tot = struct();                  % Struct per salvare tutte le FRF calcolate
field_names = fieldnames(data);    % Ottiene i nomi dei campi (dataset) in 'data'

for i = 1:length(field_names)
    Dati = data.(field_names{i}).Dati;  % Estrae la matrice dati dal file

    % Converte i dati dalla tensione alle unità fisiche
    Force_time = sens_force * Dati(:,1);
    Acc_time = sens_acc * Dati(:,2);

    % Definisce la finestra di analisi e la finestra di Hanning
    SecPoints = round(Time * fsamp);
    N_OL = floor(ol * SecPoints);
    Win = hanning(SecPoints);

    % Calcola lo spettro incrociato forza-accelerazione
    [Gxy, frq] = autocross(Force_time, Acc_time, fsamp, SecPoints, N_OL, Win);
    % Calcola lo spettro di potenza della forza
    [Gxx, ~] = autocross(Force_time, Force_time, fsamp, SecPoints, N_OL, Win);
    
    % Calcola la FRF H1
    H1 = Gxy ./ Gxx;

    % Salva la FRF calcolata nella struct
    H1_tot.(field_names{i}) = H1;
end

% Prepara le etichette per la legenda sostituendo '_' con '.'
legend_entries = strrep(field_names, '_', '.');

%% ========== PLOT PRIMO MODO (5-50 Hz) ==========

figure('Name', 'FRF - Primo Modo (5-50 Hz)', 'Position', [100 100 1200 600]);
tiledlayout(2,1);

% MODULO della FRF
nexttile;
hold off;           % Pulisce ogni plot precedente in questo axes
grid on;
title('|FRF| (5-50 Hz)');
ylabel('|H1| [m/(s^2*N)]');
xlabel('Frequenza [Hz]');

for i = 1:length(field_names)
    valid_idx = frq > f_start_plot & frq < 50;
    plot(frq(valid_idx), abs(H1_tot.(field_names{i})(valid_idx)));
    hold on;         % Mantiene il plot per aggiungere le curve successive
end
legend(legend_entries, 'Location', 'best');

% FASE della FRF
nexttile;
hold off;           % Importante: resetta la figura per evitare sovrapposizioni!
grid on;
title('Fase FRF (5-50 Hz)');
ylabel('Fase [rad]');
xlabel('Frequenza [Hz]');

for i = 1:length(field_names)
    valid_idx = frq > f_start_plot & frq < 50;
    plot(frq(valid_idx), angle(H1_tot.(field_names{i})(valid_idx)));
    hold on;
end

%% ========== PLOT SECONDO MODO (200-400 Hz) ==========

figure('Name', 'FRF - Secondo Modo (200-400 Hz)', 'Position', [150 150 1200 600]);
tiledlayout(2,1);

% MODULO della FRF
nexttile;
hold off;
grid on;
title('|FRF| (200-400 Hz)');
ylabel('|H1| [m/(s^2*N)]');
xlabel('Frequenza [Hz]');

for i = 1:length(field_names)
    valid_idx = frq > 200 & frq < 400;
    plot(frq(valid_idx), abs(H1_tot.(field_names{i})(valid_idx)));
    hold on;
end
legend(legend_entries, 'Location', 'best');

% FASE della FRF
nexttile;
hold off;
grid on;
title('Fase FRF (200-400 Hz)');
ylabel('Fase [rad]');
xlabel('Frequenza [Hz]');

for i = 1:length(field_names)
    valid_idx = frq > 200 & frq < 400;
    plot(frq(valid_idx), angle(H1_tot.(field_names{i})(valid_idx)));
    hold on;
end

%% ========== PULIZIA WORKSPACE ==========

clearvars -except data frq H1_tot
