function accIndex = selectAccelerometer(num_acc, data)
%SELEZIONAACCELEROMETRO Mostra una GUI per scegliere un accelerometro tra 1 e num_acc
%   Restituisce l'indice selezionato (intero) oppure [] se la finestra viene chiusa senza selezione.

    accIndex = [];

    % === CREA FIGURA ===
    fig = uifigure('Name', 'Seleziona Accelerometro', 'Position', [500 400 500 200]);
    fig.CloseRequestFcn = @(src, event) delete(fig);  % Permette chiusura pulita

    % === TESTO INFORMATIVO ===
    uilabel(fig, ...
        'Text', sprintf('Numero massimo accelerometri disponibili: %d', num_acc), ...
        'Position', [100 140 300 30], ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');

    % === DROPDOWN ===
    dd = uidropdown(fig, ...
        'Items', arrayfun(@(i) sprintf('Accelerometro %d', i), 1:num_acc, 'UniformOutput', false), ...
        'Position', [100 90 300 30]);

    % === BOTTONE CONFERMA ===
    uibutton(fig, ...
        'Text', 'Conferma', ...
        'Position', [200 40 100 30], ...
        'ButtonPushedFcn', @(btn, ~) confermaCallback());

    % === BLOCCA E ASPETTA SCELTA ===
    uiwait(fig);

    % === CALLBACK ===
    function confermaCallback()
        % Ottieni indice dell'accelerometro selezionato
        selectedIndex = dd.Value;
        accIndex = sscanf(selectedIndex, 'Accelerometro %d');

        % Chiudi la finestra di selezione accelerometro
        close(fig);

        % Estrai FRF dell'accelerometro selezionato
        FRF = data.frf(:, accIndex);
        f = data.fvec;

        % Salva FRF e f in base workspace per accesso nel resto dello script
        assignin('base', 'FRF', FRF);
        assignin('base', 'f', f);


        % Avvia selezione range (emulando la logica successiva già presente)
        choice = questdlg('Vuoi selezionare un range specifico della FRF?', ...
            'Selezione range', 'Sì','No','No');

        if strcmp(choice, 'Sì')
            figure('Name', 'Seleziona Range FRF', 'Position', [100 100 1600 900]);
            semilogy(f, abs(FRF)); grid on; grid minor;
            title('Seleziona due punti per il range');
            [xsel, ~] = ginput(2); close;
            fmin = min(xsel); fmax = max(xsel);
        else
            fmax = max(f);
            if min(f) == 0
                fmin = realmin;
            else
                fmin = min(f);
            end
        end

        % Filtra i dati per il range selezionato
        idx = f >= fmin & f <= fmax;
        f_range = f(idx);
        FRF_range = FRF(idx);

        % Salva anche questi nel workspace base per uso successivo
        assignin('base', 'f_range', f_range);
        assignin('base', 'FRF_range', FRF_range);

    end
end
