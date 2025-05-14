function G = TransferFunctionComputation(data)
    %TRANSFERFUNCTIONCOMPUTATION Calcola le funzioni di trasferimento (FRF),
    % le frequenze naturali e le forme modali di una trave cilindrica cava
    % a sbalzo con massa concentrata in punta, tramite metodo modale.
    
    %% === Parametri geometrici e materiali ===
    density  = data.theta_m;                % Densità materiale [kg/m^3]
    L     = data.L;                         % Lunghezza della trave [m]
    E     = data.E * 1e9;                   % Modulo di Young [Pa]
    D_ext = data.D_ext;                     % Diametro esterno [m]
    D_int = data.D_int;                     % Diametro interno [m]
    M_a   = data.M;                         % Massa concentrata in punta [kg]
    xsi = data.xi1;                         % Smorzamento critico


    % Area e momento d'inerzia della sezione trasversale
    A = (pi / 4) * (D_ext^2 - D_int^2);         
    J = (pi / 64) * (D_ext^4 - D_int^4);   

    %% === Frequenze di calcolo ===
    f_max = 700;               % Frequenza massima [Hz]
    load("time.mat"); 
    n_points = length(t); % Numero di punti 
    freq = linspace(0, f_max, n_points);     % Vettore di frequenze in Hz (0 -> f_max)
    omega = 2 * pi * freq;            % Vettore di pulsazioni in rad/s


    %% === Calcolo del determinante della matrice e frequenze naturali ===
    dets = MatrixDeterminant(omega, density, A, E, J, M_a, L);        
    [i_nat, G] = NaturalFrequencyComputation(dets, freq);       
    
    %% === Calcolo coefficienti modali ===
    X_hat = SystemSolver(i_nat, omega, density, A, E, J, M_a, L);
    
    %% === Costruzione delle forme modali ===
    [x, modes_shapes] = ModeShapes(L, n_points, i_nat, omega, density, A, E, J, X_hat);
    
    %% === Definizione posizione forza e sensore ===
    xj = data.L;  % Posizione accelerometro [m]

    [~, pos_xj] = min(abs(x - xj));  
    
    %% === Calcolo massa modale ===
    modal_mass = zeros(1, length(i_nat));
    for i = 1:length(i_nat)
        phi = modes_shapes(i,:);  % forma modale
        % massa distribuita integrata
        int_mass = trapz(x, density * A * phi.^2);  
        % massa concentrata alla fine
        tip_mass = M_a * phi(end)^2;  
        % massa modale totale
        modal_mass(i) = int_mass + tip_mass;
    end

    %% === Calcolo delle FRF: G_WY e G_SY ===
    [phi_at_xj, G] = GwyComputation(omega, i_nat, modes_shapes, x, pos_xj, modal_mass, xsi, density, A, M_a, G);
    G = GsyComputation(omega, i_nat, modes_shapes, modal_mass, xsi, phi_at_xj, G, E, J, x, pos_xj, density, A, M_a);

    %% === SALVO I VETTORI DI FREQ E OMEGA ===
    G.freq = freq;
    G.omega = omega;
    
    %% === Plot finale ===
    GPlotSingle(freq, G);  % Visualizzazione risultati
end





