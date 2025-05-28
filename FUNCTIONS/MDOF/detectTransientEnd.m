function kEnd = detectTransientEnd(IRF, fracThreshold, windowSize)
    % IRF: Nt x nChan x nAcq
    % time: Nt x 1
    % fracThreshold: soglia relativa (es. 0.05 = 5%)
    % windowSize: numero di campioni consecutivi per validare il taglio

    [Nt, ~, ~] = size(IRF.irf);
    time = IRF.time;

    % 1) Calcolo ampiezza RMS globale A(k)
    A = zeros(Nt,1);
    for k = 1:Nt
        block = IRF.irf(k,:,:);
        A(k) = sqrt( mean( abs(block(:)).^2 ) );
    end

    % 2) Soglia assoluta
    Amax = max(A);
    Ath = fracThreshold * Amax;

    % 3) Trova primo k dove A(k) < Ath per windowSize campioni
    below = A < Ath;
    kEnd = Nt;  % default: tutta la lunghezza
    for k = 1:Nt-windowSize+1
        if all(below(k:k+windowSize-1))
            kEnd = k;
            break;
        end
    end

    % Visualizzazione (opzionale)
    figure;
    plot(time, A, '-'); hold on;
    yline(Ath, '--r', sprintf('%.1f%% A_{max}', fracThreshold*100));
    xline(time(kEnd), '--k', 'Cut-off');
    xlabel('Time (s)'); ylabel('Global RMS amplitude');
    title('Rilevamento fine primo transitorio');
    grid on;
end
