function [stable_freqs] = stabilization_diagram_old(IRF, n_range, kEnd,tol_freq, tol_damp)
% build_stabilization_diagram Costruisce il diagramma di stabilizzazione
%
%   stabilization_diagram(IRF, time, n_min, n_max, m, tol_freq, tol_damp)
%   IRF: struct con IRF.irf_cut (Nt x nChan x nAcq)
%   time: vettore tempo
%   n_min, n_max: intervallo ordini
%   m: numero righe Hankel
%   tol_freq: tolleranza frequenza (es. 0.01 → 1%)
%   tol_damp: tolleranza damping (es. 0.05 → 5%)
    time = IRF.time;
    dt = time(2) - time(1);
    stable_freqs = [];

    figure; hold on;
    title('Stabilization Diagram');
    xlabel('Order n');
    ylabel('Frequency [Hz]');

    for n = n_range
        m = kEnd - n;
        [H0, H1] = build_hankel(IRF, m, n);
        H0_pinv = compute_pseudoinverse(H0);
        A = H1 * H0_pinv;
        [freq, damping] = extract_modal_parameters(A, dt);

        % Plot tutti i modi (grigi)
        plot(n * ones(size(freq)), freq, 'ko', 'MarkerSize', 3);

        % Verifica stabilità rispetto al passo precedente
        if n > n_min
            for i = 1:length(freq)
                f_curr = freq(i);
                d_curr = damping(i);

                for j = 1:size(stable_freqs,1)
                    f_prev = stable_freqs(j,1);
                    d_prev = stable_freqs(j,2);

                    if abs(f_curr - f_prev)/f_prev < tol_freq && ...
                       abs(d_curr - d_prev)/d_prev < tol_damp
                        % Stabile
                        plot(n, f_curr, 'bo', 'MarkerFaceColor', 'b');
                        break;
                    end
                end
            end
        end

        % Salva modi correnti
        stable_freqs = [freq(:), damping(:)];
    end

    grid on;
end
