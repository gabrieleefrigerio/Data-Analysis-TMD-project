function [data] = PowerBalance(data)
    
     % === Descrizione della funzione ===
    % La funzione PowerBalance calcola il bilancio di potenza per una trave
    % di materiale SMA (Shape Memory Alloy) sottoposta a corrente elettrica.
    % La potenza generata dall'effetto Joule (P_joule) viene confrontata con
    % la potenza dissipata per convezione (P_conv), e il sistema raggiunge una
    % temperatura di equilibrio (T_eq). La funzione aggiorna la temperatura
    % della trave in base a questi bilanci energetici, e restituisce i risultati
    % come la resistenza elettrica (R), la potenza generata per effetto Joule
    % (P_joule), la potenza dissipata per convezione (P_conv), e la differenza
    % tra le due potenze (P_diff), che deve tendere a zero per un sistema in equilibrio.
    %
    % === Parametri di input ===
    % data: una struttura contenente i dati relativi alla trave, inclusi:
    % - Corrente (Current)
    % - Resistenza del materiale (rho)
    % - Coefficiente di scambio termico (h)
    % - Temperature e proprietà dell'aria
    %
    % === Parametri di output ===
    % data: la struttura aggiornata contenente i risultati del bilancio energetico,
    % tra cui la temperatura aggiornata della trave (T), la resistenza (R),
    % e la potenza dissipata (P_joule, P_conv, P_diff).

    % === Corrente elettrica nella trave ===
    I = data.Current;  % [A]

    % === Parametri geometrici ===
    D_ext = data.D_ext;      % [m] Diametro esterno
    D_int = data.D_int;      % [m] Diametro interno
    L     = data.L;          % [m] Lunghezza della trave

    % === Area della sezione trasversale cava ===
    A_sec = pi/4 * (D_ext^2 - D_int^2);  % [m^2]

    % === Superficie esposta alla convezione (esterna + interna) ===
    A_surf = pi * (D_ext + D_int) * L;   % [m^2]

    % === Calcolo resistenza elettrica ===
    R = data.rho * L / A_sec;           % [Ohm]

    % === Potenza generata per effetto Joule ===
    P_joule = I^2 * R;                  % [W]

    % === Calcolo temperatura a regime ===
    % P_joule = h * A_surf * (T - T_aria) --> T = P_joule / (h * A_surf) + T_aria
    h = data.h;                         % [W/(m²·K)]
    T_air = data.air.Tair;              % [K]
    T_eq = P_joule / (h * A_surf) + T_air;  % [K]

    % === Salvo la temperatura ===
    data.Tprev = data.T;
    data.T     = T_eq;

    % === Potenza dissipata (coerente con la T trovata) ===
    P_conv = h * A_surf * (T_eq - T_air);  % [W]

    % === Salvataggio risultati ===
    data.R        = R;
    data.P_joule  = P_joule;
    data.P_conv   = P_conv;
    data.P_diff   = P_joule - P_conv;  % Deve tendere a zero

end
