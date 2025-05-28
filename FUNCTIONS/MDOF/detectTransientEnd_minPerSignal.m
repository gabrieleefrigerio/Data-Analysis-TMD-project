function kEnd = detectTransientEnd_minPerSignal(IRF, fracThreshold, windowSize)
    % IRF: Nt x nChan x nAcq
    % fracThreshold: soglia relativa (es. 0.05 = 5%)
    % windowSize: numero di campioni consecutivi per dire che è finito

    [Nt, nChan, nAcq] = size(IRF.irf);
    time = IRF.time;
    Amax_all = squeeze(max(abs(IRF.irf), [], 1));  % max per ogni (chan, acq)
     if nAcq == 1
        Amax_all = Amax_all';
    end
    kEnd_all = Nt * ones(nChan, nAcq);  % inizializza a valore massimo

    for i = 1:nChan
        for j = 1:nAcq
            signal = squeeze(IRF.irf(:,i,j));
            A = abs(signal);  % oppure RMS con filtro, se vuoi

            Amax = Amax_all(i,j);
            Ath = fracThreshold * Amax;

            below = A < Ath;

            for k = 1:Nt - windowSize + 1
                if all(below(k:k + windowSize - 1))
                    kEnd_all(i,j) = k;
                    break;
                end
            end
        end
    end

    kEnd = min(kEnd_all(:));  % il più piccolo kEnd tra tutti

    % Visualizza
    figure; hold on;
    plot(time, squeeze(mean(mean(abs(IRF.irf), 2), 3)), 'b');
    xline(time(kEnd), '--r', 'Global min k_{end}');
    xlabel('Time (s)');
    ylabel('Average amplitude');
    title('Rilevamento k_{end} per ogni segnale (min globale)');
    grid on;
end
