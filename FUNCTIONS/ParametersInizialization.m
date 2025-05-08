function data = ParametersInizialization()
    % === Descrizione della funzione ===
    % La funzione ParametersInizialization ha il compito di inizializzare
    % i parametri fisici e geometrici necessari per il modello di una trave
    % in materiale SMA (Shape Memory Alloy) utilizzata in un TMD (Tuned Mass Damper).
    % I parametri sono forniti in base alle proprietà specifiche del materiale (NiTiNOL)
    % e alle condizioni ambientali come la temperatura dell'aria. La funzione
    % restituisce una struttura 'data' con tutti i parametri inizializzati.
    %
    % === Parametri di output ===
    % data: la struttura aggiornata con i parametri di configurazione necessari
    % per il modello SMA e per i calcoli successivi (come la temperatura, la
    % resistività, le proprietà termiche e geometriche della trave).

    % === Inizializzazione della struttura 'data' ===
    data = struct();

    % === PAPER PARAMETERS ===
    % Parametri del primo materiale - SMA TMD
    % Inizializzazione dei parametri che definiscono le caratteristiche del
    % materiale SMA (Shape Memory Alloy), inclusi i punti di transizione tra
    % Austenite e Martensite, i moduli elastici e la resistività.
    data.As      = 55.0;        % [°C] Austenite start
    data.Af      = 65.0;        % [°C] Austenite finish
    data.Ms      = 40.0;        % [°C] Martensite start
    data.Mf      = 28.5;        % [°C] Martensite finish

    data.Em      = 32.3;        % [GPa] Modulo elastico Martensite
    data.Ea      = 52.7;        % [GPa] Modulo elastico Austenite

    data.rho_m   = 90e-8;       % [Ohm·m] Resistività Martensite
    data.rho_a   = 100e-8;      % [Ohm·m] Resistività Austenite

    data.xi1_m = 1.22e-2;       % [-] Smorzamento Martensite
    data.xi1_a = 0.90e-2;       % [-] Smorzamento Austenite

    data.T = 25;                % [°C] Temperatura Tmd iniziale
    data.T_prev = data.T;       % [°C] Temperatura Tmd precedente

    % === Proprietà dell'aria (a 25°C, 1 atm) ===
    % Inizializzazione delle proprietà dell'aria a temperatura ambiente.
    % Queste proprietà sono necessarie per calcolare il raffreddamento per convezione
    % della trave SMA in base alla temperatura dell'aria e al coefficiente di scambio termico.
    data.air.Tair = 25;         % [°C] Temperatura aria
    data.air.rho_a   = 1.184;   % [kg/m^3] Densità dell'aria
    data.air.ba      = 1.849e-5; % [Pa·s] Viscosità dinamica dell'aria
    data.air.eta_a   = 1005;    % [J/(kg·K)] Calore specifico a pressione costante
    data.air.fi_a    = 0.0262;  % [W/(m·K)] Conducibilità termica

    % === Proprietà geometriche della trave SMA (NiTiNOL) ===
    % Inizializzazione delle proprietà geometriche della trave SMA,
    % inclusi la lunghezza, i diametri e la massa concentrata.
    data.L      = 0.14;         % [m] Lunghezza della trave (140 mm)
    data.D_ext  = 0.004;        % [m] Diametro esterno della sezione circolare cava (4 mm)
    data.D_int  = 0.003;        % [m] Diametro interno, dato da D_ext - 2*t (3 mm, con spessore 0.5 mm)
    data.M      = 0.0176;       % [kg] Massa concentrata applicata a ciascuna estremità

    % === Parametri opzionali e placeholder ===
    % Alcuni parametri, come i coefficienti termomeccanici e l'espansione termica,
    % non sono forniti nel materiale di riferimento e sono lasciati come placeholder.
    data.CA      = NaN;         % [MPa/°C] Coefficiente termomeccanico Austenite
    data.CM      = NaN;         % [MPa/°C] Coefficiente termomeccanico Martensite
    data.Hcur    = NaN;         % [-] Parametro curva isteresi
    data.alpha   = NaN;         % [°C^-1] Coeff. di espansione termica

end

