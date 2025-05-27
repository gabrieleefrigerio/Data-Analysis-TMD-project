function GPlot(freq, I, GSY_values, GWY_values)
    % GPlot: Plotta ampiezza e fase di G_SY e G_WY su 4 subplot verticali con asse Y logaritmico per l'ampiezza.

    % Crea figura docked
    fig = figure('Name', 'FRF SMA TMD', 'NumberTitle', 'off', 'WindowStyle', 'docked');
    
    % Ampiezza G_SY
    subplot(4, 1, 1);
    hold on;
    for ii = 1:length(I)
        plot(freq, abs(GSY_values(:,ii)), 'LineWidth', 2);
    end
    set(gca, 'YScale', 'log');
    ylabel('|G_{SY}|');
    title('Ampiezza G_{SY} vs Frequenza');
    legend(arrayfun(@(x) sprintf('I = %.1f A', x), I, 'UniformOutput', false), 'Location', 'northeastoutside');
    grid on;
    hold off;

    % Fase G_SY
    subplot(4, 1, 2);
    hold on;
    for ii = 1:length(I)
        plot(freq, unwrap(angle(GSY_values(:,ii))) * 180/pi, 'LineWidth', 2);
    end
    ylabel('Fase G_{SY} [°]');
    title('Fase G_{SY} vs Frequenza');
    legend('off');
    grid on;
    hold off;

    % Ampiezza G_WY
    subplot(4, 1, 3);
    hold on;
    for ii = 1:length(I)
        plot(freq, abs(GWY_values(:,ii)), 'LineWidth', 2);
    end
    set(gca, 'YScale', 'log');
    ylabel('|G_{WY}|');
    title('Ampiezza G_{WY} vs Frequenza');
    legend('off');
    grid on;
    hold off;

    % Fase G_WY
    subplot(4, 1, 4);
    hold on;
    for ii = 1:length(I)
        plot(freq, unwrap(angle(GWY_values(:,ii))) * 180/pi, 'LineWidth', 2);
    end
    xlabel('Frequenza [Hz]');
    ylabel('Fase G_{WY} [°]');
    title('Fase G_{WY} vs Frequenza');
    legend('off');
    grid on;
    hold off;

    % Ottimizza aspetto finale
    set(fig, 'Color', 'w');
end
