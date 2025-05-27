function fitSingleMode(modeIndex, peaks, locs, f_range, FRF_range, allParams, acqNumber, accIndex)
    if modeIndex > length(peaks)
        showFinalResult(allParams, f_range, FRF_range, accIndex, acqNumber);
        return;
    end

    figFit = uifigure('Name', ['Modo ', num2str(modeIndex)], 'Position', [100 100 900 600]);

    f0 = locs(modeIndex);
    [~, locs_idx] = min(abs(f_range - f0));

    mag_target = abs(FRF_range(locs_idx))/sqrt(2);
    left_idx = find(abs(FRF_range(1:locs_idx)) <= mag_target, 1, 'last');
    right_idx = find(abs(FRF_range(locs_idx:end)) <= mag_target, 1, 'first') + locs_idx - 1;

    if isempty(left_idx) || isempty(right_idx)
        xi0 = 0.01;
    else
        f1 = f_range(left_idx);
        f2 = f_range(right_idx);
        xi0 = (f2 - f1) / (2*f0);
    end

    omega0 = 2*pi*f0;
    A0 = real(FRF_range(locs_idx) * (2i * xi0 * omega0^2));
    Rh0 = 0;
    Rl0 = 0;

    fieldf0 = uieditfield(figFit, 'numeric', 'Value', f0, 'Position', [150 450 100 22]);
    uilabel(figFit, 'Text', 'f0 (Hz):', 'Position', [30 450 100 22]);

    fieldXi = uieditfield(figFit, 'numeric', 'Value', xi0, 'Position', [150 410 100 22]);
    uilabel(figFit, 'Text', 'xi:', 'Position', [30 410 100 22]);

    fieldA = uieditfield(figFit, 'numeric', 'Value', A0, 'Position', [150 370 100 22]);
    uilabel(figFit, 'Text', 'A:', 'Position', [30 370 100 22]);

    fieldRh = uieditfield(figFit, 'numeric', 'Value', Rh0, 'Position', [150 330 100 22]);
    uilabel(figFit, 'Text', 'Rh:', 'Position', [30 330 100 22]);

    fieldRl = uieditfield(figFit, 'numeric', 'Value', Rl0, 'Position', [150 290 100 22]);
    uilabel(figFit, 'Text', 'Rl:', 'Position', [30 290 100 22]);

    btnOptimize = uibutton(figFit, 'Text', 'Ottimizza', ...
            'Position', [30 200 220 30], ...
            'ButtonPushedFcn', @(~,~) runOptimization());

    btnNext = uibutton(figFit, 'Text', 'Avanti', ...
            'Position', [30 160 220 30], ...
            'ButtonPushedFcn', @(~,~) proceedToNext());

    axAmp = uiaxes(figFit, 'Position', [300, 350, 580, 250]);
    title(axAmp, 'FRF - Ampiezza');
    xlim(axAmp, [f0 - 50, f0 + 50]);
    xlabel(axAmp, 'Frequenza (Hz)');
    ylabel(axAmp, '|FRF|');
    grid(axAmp, 'on'); grid(axAmp, 'minor');

    axPhase = uiaxes(figFit, 'Position', [300, 50, 580, 250]);
    title(axPhase, 'FRF - Fase');
    xlim(axPhase, [f0 - 50, f0 + 50]);
    xlabel(axPhase, 'Frequenza (Hz)');
    ylabel(axPhase, 'Fase [rad]');
    grid(axPhase, 'on'); grid(axPhase, 'minor');
    optimized = [];
  runOptimization();

    %% FUNZIONE OTTIMIZZAZIONE
    function runOptimization()
            % vettore dei parametri da ottimizzare
            p0 = [fieldf0.Value*2*pi, fieldXi.Value, fieldA.Value, fieldRh.Value, fieldRl.Value];

            % === Selezione dinamica dell'intervallo attorno al picco ===

            % Frequenza centrale del picco corrente
            f_central = locs(modeIndex);

            % Calcolo della larghezza della finestra: metà distanza tra picchi adiacenti
            resol = 15;
            if modeIndex == 1
                if length(modeIndex) == 1
                    df = 50;
                else
                % Primo picco: guarda solo verso il prossimo
                df = (locs(2) - locs(1)) / resol;
                end
            elseif modeIndex == length(locs)
                % Ultimo picco: guarda solo verso il precedente
                df = (locs(end) - locs(end-1)) / resol;
            else
                % Picchi centrali: usa la media tra le due distanze adiacenti
                df1 = locs(modeIndex) - locs(modeIndex - 1);
                df2 = locs(modeIndex + 1) - locs(modeIndex);
                df = min(df1, df2) / resol; % più conservativo
            end

            % Estendi di un piccolo fattore (es. 20%) per sicurezza
            f_low = f_central -  df;
            f_high = f_central +  df;

            % Trova gli indici nell'intervallo
            idx_min = find(f_range >= f_low, 1, 'first');
            idx_max = find(f_range <= f_high, 1, 'last');

            % Fallback agli estremi se gli indici non vengono trovati
            if isempty(idx_min), idx_min = 1; end
            if isempty(idx_max), idx_max = length(f_range); end

            % Vettori finali da usare
            omega_vec = 2 * pi * f_range(idx_min:idx_max);
            G_exp = FRF_range(idx_min:idx_max).';


            % funzione di trasferimento numerica in formato anonymous
            modelFun = @(p, omega_vec) p(3)./ (-omega_vec.^2 + 2j*p(2)*p(1).*omega_vec + p(1)^2) + p(4) + ( p(5)./omega_vec.^2);
            % funzione con parametri scalati
            scale = [1, 1, 1, 1e-3, 1e-6];
            modelFun_scaled = @(p, omega_vec) modelFun(p .* scale, omega_vec).';

            % cost function da minimizzare
            residui = @(p) sum( real(G_exp - modelFun_scaled(p, omega_vec)).^2 + imag(G_exp - modelFun_scaled(p, omega_vec)).^2  );

            % setting lsqnonlin
            opts = optimoptions('lsqnonlin', 'Display', 'iter', 'TolFun',1e-12, 'TolX',1e-12);

            % effettuo l'ottimizzazione
            [popt_scaled, ~] = lsqnonlin(residui, p0, [], [], opts);

            popt = popt_scaled .* scale;

            % salvo i valori ottimizzati
            optimized = popt;

            % calcolo la funzione di trasf numerica con i parmaetri
            % ottimizzati
            G_fit = modelFun(popt, omega_vec);

            % Aggiorno il plot
            % ---- Ampiezza ----
            cla(axAmp);
            semilogy(axAmp, f_range(idx_min:idx_max), abs(G_exp), 'b'); hold(axAmp, 'on'); grid(axAmp, 'minor');% grid(axAmp, 'on');
            xlim(axAmp,[fieldf0.Value-50 fieldf0.Value+50]) % plotto il grafico solo nell'intorno del picco
            semilogy(axAmp, f_range(idx_min:idx_max), abs(G_fit), 'r--', 'LineWidth', 1.5);
            legend(axAmp, 'FRF', 'Fit');
            grid(axAmp, 'minor'); %grid(axAmp, 'on');
            hold(axAmp, 'off');

            % ---- Fase ----
            cla(axPhase);
            plot(axPhase, f_range(idx_min:idx_max), angle(G_exp), 'b'); hold(axPhase, 'on');  grid(axPhase, 'minor'); % grid(axPhase, 'on');
            xlim(axPhase,[fieldf0.Value-50 fieldf0.Value+50]) % plotto il grafico solo nell'intorno del picco
            plot(axPhase, f_range(idx_min:idx_max), angle(G_fit), 'r--', 'LineWidth', 1.5);
            legend(axPhase, 'FRF', 'Fit');
            grid(axPhase, 'minor'); %grid(axPhase, 'on');
            hold(axPhase, 'off');

            % Aggiorna i campi nella gui
            fieldf0.Value = popt(1)/2/pi;
            fieldXi.Value = popt(2);
            fieldA.Value = popt(3);
            fieldRh.Value = popt(4);
            fieldRl.Value = popt(5);
        end
 %% FUNZIONE TASTO PER PROCEDERE AL OTTIMIZZAZIONE DOPO
        function proceedToNext()
            if isempty(optimized)
                runOptimization();
            end
            % chiudo questa figure
            close(figFit);

            % apro la figure per il modo seguente
            fitSingleMode(modeIndex + 1, peaks, locs, f_range, FRF_range, [allParams; optimized], acqNumber, accIndex);
        end
end
