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
data.I5_8 = load('Corrente_5.8A\Prova2_4min\test@3830842267.mat');
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

H1_tot = struct();                  % Struct per le FRF
coherence_tot = struct();          % Struct per la coerenza
field_names = fieldnames(data);

for i = 1:length(field_names)
    Dati = data.(field_names{i}).Dati;

    Force_time = sens_force * Dati(:,1);
    Acc_time = sens_acc * Dati(:,2);

    SecPoints = round(Time * fsamp);
    N_OL = floor(ol * SecPoints);
    Win = hanning(SecPoints);
    [Gxy, frq] = autocross_lab(Acc_time, Force_time, fsamp, SecPoints, N_OL, Win);
    [Gyy, ~] = autocross_lab(Force_time, Force_time, fsamp, SecPoints, N_OL, Win);
    [Gxx, ~] = autocross_lab(Acc_time, Acc_time, fsamp, SecPoints, N_OL, Win);
    
    omega = (frq/(2*pi))';
    %H1 = ((Gxy+0.88).*(-(omega).^2) ./ Gxx);


    H1 = (Gxy./ Gxx)+0.88;
    H2 = (Gxy) ./ Gyy;
    %gamma2 = abs(Gxy).^2 ./ (Gxx .* Gyy);
    gamma2 = abs(Gxy).^2 ./ (Gxx .* Gyy);

    H1_tot.(field_names{i}) = H1;
    coherence_tot.(field_names{i}) = gamma2;
end


% Prepara le etichette per la legenda sostituendo '_' con '.'
legend_entries = strrep(field_names, '_', '.');

% Bande di frequenza per il primo e secondo modo
f_start_plot_1 = 5;
f_finish_plot_1 = 50;
f_start_plot_2 = 200;
f_finish_plot_2 = 400;

%% ========== PLOT PRIMO MODO (5-50 Hz) ==========

figure('Name', 'FRF - Primo Modo (5-50 Hz)', 'Position', [50 50 1200 800]);
tiledlayout(3,1);

% Modulo |H1|
nexttile;
hold on; grid on; grid minor;
title('Modulo |FRF| - Primo Modo (5-50 Hz)');
xlabel('Frequenza [Hz]');
ylabel('|H1| [N/m]');
for i = 1:length(field_names)
    valid_idx = frq > f_start_plot_1 & frq < f_finish_plot_1;
    plot(frq(valid_idx), abs(H1_tot.(field_names{i})(valid_idx)),LineWidth=1.2);
end
legend(legend_entries, 'Location', 'best');

% Fase
nexttile;
hold on; grid on; grid minor;
title('Fase FRF - Primo Modo (5-50 Hz)');
xlabel('Frequenza [Hz]');
ylabel('Fase [rad]');
for i = 1:length(field_names)
    valid_idx = frq > f_start_plot_1 & frq < f_finish_plot_1;
    plot(frq(valid_idx), angle(H1_tot.(field_names{i})(valid_idx)),LineWidth=1.5);
end

% Coerenza
nexttile;
hold on; grid on; grid minor;
title('Coerenza - Primo Modo (5-50 Hz)');
xlabel('Frequenza [Hz]');
ylabel('\gamma^2');
ylim([0 1.05]);
for i = 1:length(field_names)
    valid_idx = frq > f_start_plot_1 & frq < f_finish_plot_1;
    plot(frq(valid_idx), coherence_tot.(field_names{i})(valid_idx), LineWidth=1.2);
end

%% ========== PLOT SECONDO MODO (200-400 Hz) ==========

figure('Name', 'FRF - Secondo Modo (200-400 Hz)', 'Position', [50 50 1200 800]);
tiledlayout(3,1);

% Modulo
nexttile;
hold on; grid on; grid minor;
title('Modulo |FRF| - Secondo Modo (200-400 Hz)');
xlabel('Frequenza [Hz]');
ylabel('|H1| [N/m]');
for i = 1:length(field_names)
    valid_idx = frq > f_start_plot_2 & frq < f_finish_plot_2;
    plot(frq(valid_idx), abs(H1_tot.(field_names{i})(valid_idx)),LineWidth=1.5);
end
legend(legend_entries, 'Location', 'best');

% Fase
nexttile;
hold on; grid on; grid minor;
title('Fase FRF - Secondo Modo (200-400 Hz)');
xlabel('Frequenza [Hz]');
ylabel('Fase [rad]');
for i = 1:length(field_names)
    valid_idx = frq > f_start_plot_2 & frq < f_finish_plot_2;
    plot(frq(valid_idx), angle(H1_tot.(field_names{i})(valid_idx)));
end

% Coerenza
nexttile;
hold on; grid on; grid minor;
title('Coerenza - Secondo Modo (200-400 Hz)');
xlabel('Frequenza [Hz]');
ylabel('\gamma^2');
ylim([0 1.05]);
for i = 1:length(field_names)
    valid_idx = frq > f_start_plot_2 & frq < f_finish_plot_2;
    plot(frq(valid_idx), coherence_tot.(field_names{i})(valid_idx));
end

%% ========== PLOT COMPLETO LOG-SCALE ==========

figure('Name', 'FRF Completa - Scala logaritmica', 'Position', [50 50 1200 800]);
tiledlayout(3,1);

% Modulo in scala logaritmica
nexttile;
hold on; grid on; grid minor;
title('Modulo |FRF| - Banda Completa');
xlabel('Frequenza [Hz]');
ylabel('|H1| [N/m]');
set(gca, 'YScale', 'log');
for i = 1:length(field_names)
    valid_idx = frq > f_start_plot_1 & frq < f_finish_plot_2;
    semilogy(frq(valid_idx), abs(H1_tot.(field_names{i})(valid_idx)),LineWidth=1.5);
end
legend(legend_entries, 'Location', 'best');

% Fase
nexttile;
hold on; grid on; grid minor;
title('Fase FRF - Banda Completa');
xlabel('Frequenza [Hz]');
ylabel('Fase [rad]');
for i = 1:length(field_names)
    valid_idx = frq > f_start_plot_1 & frq < f_finish_plot_2;
    plot(frq(valid_idx), angle(H1_tot.(field_names{i})(valid_idx)));
end

% Coerenza
nexttile;
hold on; grid on; grid minor;
title('Coerenza - Banda Completa');
xlabel('Frequenza [Hz]');
ylabel('\gamma^2');
ylim([0 1.05]);
for i = 1:length(field_names)
    valid_idx = frq > f_start_plot_1 & frq < f_finish_plot_2;
    plot(frq(valid_idx), coherence_tot.(field_names{i})(valid_idx));
end

%% ========== ESPORTAZIONE FRF (TAGLIO A 600 Hz) IN FILE .mat SEPARATI ==========

output_folder = 'exp_for_ibrahim';  % Cartella di destinazione per i file esportati

% Crea la cartella se non esiste
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Filtro frequenze <= 600 Hz
freq_cut_idx = frq <= 600;
frq = frq(freq_cut_idx);  % Sovrascrive con versione tagliata

% Itera su ciascun campo della struttura H1_tot
for i = 1:length(field_names)
    field = field_names{i};
    
    % Applica il filtro alle FRF e alla coerenza
    H1 = H1_tot.(field);
    gamma2 = coherence_tot.(field);
    
    H1 = H1(freq_cut_idx);
    gamma2 = gamma2(freq_cut_idx);

    % Prepara il nome del file sostituendo '_' con '.' per indicare la corrente
    current_label = strrep(field, '_', '.');
    filename = fullfile(output_folder, ['FRF corrente ', current_label, '.mat']);
    
    % Salva le variabili tagliate con nomi compatibili
    save(filename, 'H1', 'gamma2', 'frq');
end


%% ========== PULIZIA E SALVATAGGIO ==========

% Indici delle frequenze nel range desiderato
valid_idx = frq >= f_start_plot_1 & frq <= f_finish_plot_2;
frq_filtered = frq(valid_idx);

% Crea una nuova struct contenente solo le FRF filtrate
H1_filtered = struct();
for i = 1:length(field_names)
    H1_filtered.(field_names{i}) = H1_tot.(field_names{i})(valid_idx);
end

% Salvataggio
save('FRF_results.mat', 'H1_filtered', 'frq_filtered');

% Pulisci variabili opzionalmente
clearvars -except data frq_filtered H1_filtered

