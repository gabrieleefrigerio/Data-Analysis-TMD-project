function [x, modes_shapes] = ModeShapes(L, n_points, i_nat, omega, density, A, E, J, X_hat)
% MODESHAPES Costruisce le forme modali normalizzate lungo una trave.
% 
% INPUT:
%   L        - Lunghezza della trave
%   n_points - Numero di punti spaziali in cui valutare la forma modale
%   i_nat    - Indici delle frequenze naturali nel vettore omega
%   omega    - Vettore delle pulsazioni naturali (rad/s)
%   density  - Densità del materiale della trave (ρ)
%   A        - Area della sezione trasversale
%   E        - Modulo di Young
%   J        - Momento d'inerzia della sezione
%   X_hat    - Matrice 4xN dei coefficienti modali per ogni modo (da SystemSolver)
%
% OUTPUT:
%   x            - Vettore delle coordinate spaziali
%   modes_shapes - Matrice NxM delle forme modali normalizzate (una riga per ogni modo)

    % Discretizzazione della trave nello spazio
    x = linspace(0, L, n_points);
    
    % Preallocazione della matrice delle forme modali
    modes_shapes = zeros(length(i_nat), length(x));
    
    % Cicla su ogni modo proprio per calcolare la rispettiva forma modale
    for i_mode = 1:length(i_nat)
        w = omega(i_nat(i_mode));  % Frequenza angolare del modo attuale
        gamma = sqrt(w) * (density * A / (E * J))^(1/4);  % Calcolo di gamma

        % Costruzione della forma modale per ogni punto x
        modes_shapes(i_mode, :) = ...
            X_hat(1, i_mode) * sin(gamma * x) + ...
            X_hat(2, i_mode) * cos(gamma * x) + ...
            X_hat(3, i_mode) * sinh(gamma * x) + ...
            X_hat(4, i_mode) * cosh(gamma * x);
    end
    
    % Normalizzazione delle forme modali in modo che il massimo valore assoluto sia 1
    modes_shapes = modes_shapes ./ max(abs(modes_shapes), [], 2);

    %% === ANIMAZIONE DEL PRIMO MODO ===
    % L'animazione qui sotto visualizza l'evoluzione temporale della deformazione modale
    % (Attualmente commentata, decommenta per attivarla)
    
    mode = 2;  % Seleziona il primo modo da animare
    freq = omega / (2*pi);  % Conversione da rad/s a Hz
    colors_p = lines(length(i_nat));  % Genera palette di colori

    % Creazione della figura e preparazione al plot
    figure('Position', [100, 100, 1200, 500]); 
    hold on; grid on;
    title(['Animazione del modo ', num2str(mode)])
    plot(x, modes_shapes(mode,:), ':k', 'LineWidth', 2)  % Traccia forma modale fissa
    h1 = plot(x, zeros(size(x)), 'LineWidth', 2, 'color', colors_p(mode, :));  % Traccia la deformazione animata
    xlabel('Posizione [m]')
    ylabel('Deformazione (modale)')
    ylim([-1, 1] * max(abs(modes_shapes(mode,:))) * 1.1)  % Imposta i limiti dell’asse y

    % Parametri animazione
    T = 1 / freq(i_nat(mode));   % Periodo del modo selezionato
    n_cycles = 4;                % Numero di cicli da visualizzare
    n_frames = 1000;             % Numero totale di frame
    t = linspace(0, n_cycles*T, n_frames);  % Vettore temporale

    % Loop di animazione
    for k = 1:length(t)
        if ishandle(h1)
            % Aggiorna la deformata modale nel tempo
            w1 = modes_shapes(mode,:) * cos(2*pi*freq(i_nat(mode)) * t(k));
            h1.YData = w1;
            pause(T/n_frames);  % Pausa per sincronizzare il frame rate
        else
            return  % Esce se la figura è chiusa
        end
    end
end


