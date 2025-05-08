function [data] = ConcentrationComputation(data)
    % ConcentrationComputation - Calcola il rapporto di concentrazione di martensite 
    % in un materiale a memoria di forma (SMA) in base alla temperatura attuale e 
    % a quella del passo precedente. Il rapporto varia tra 1 (solo martensite) e 
    % 0 (solo austenite), con interpolazione lineare nelle regioni di transizione.
    %
    % Il comportamento tiene conto della direzione del ciclo termico:
    %   - Se T aumenta → riscaldamento (Martensite → Austenite)
    %   - Se T diminuisce → raffreddamento (Austenite → Martensite)
    %   - Se T è costante → il materiale è a regime, la concentrazione resta invariata.
    %
    % Richiede i seguenti campi nella struttura 'data':
    %   data.T        : temperatura attuale [°C]
    %   data.T_prev   : temperatura al passo precedente [°C]
    %   data.As, Af   : temperature di inizio/fine trasformazione austenitica
    %   data.Ms, Mf   : temperature di inizio/fine trasformazione martensitica
    %
    % Restituisce:
    %   data.ratio    : concentrazione di martensite (1 = 100% martensite)
    
    if data.T > data.T_prev
        % --- Riscaldamento ---
        % Se la temperatura aumenta, il materiale si trasforma da Martensite ad Austenite
        if data.T < data.As
            data.ratio = 1;  % Solo Martensite
        elseif data.T < data.Af
            % Interpolazione lineare tra Martensite e Austenite nella zona di transizione
            data.ratio = (data.Af - data.T) / (data.Af - data.As);
        else
            data.ratio = 0;  % Solo Austenite
        end

    elseif data.T < data.T_prev
        % --- Raffreddamento ---
        % Se la temperatura diminuisce, il materiale si trasforma da Austenite a Martensite
        if data.T < data.Mf
            data.ratio = 1;  % Solo Martensite
        elseif data.T < data.Ms
            % Interpolazione lineare tra Austenite e Martensite nella zona di transizione
            data.ratio = (data.Mf - data.T) / (data.Mf - data.Ms);
        else
            data.ratio = 0;  % Solo Austenite
        end

    else
        % --- Temperatura costante: corpo a regime termico ---
        % Se la temperatura non cambia, il materiale è in stato di equilibrio termico
        % Se è il primo passo e 'ratio' non esiste, stimiamo in base alla temperatura
        if ~isfield(data, 'ratio')
            if data.T <= data.Mf
                data.ratio = 1;  % Martensite completo
            elseif data.T >= data.Af
                data.ratio = 0;  % Austenite completo
            else
                data.ratio = 0.5;  % Default nella zona di transizione (50% Martensite, 50% Austenite)
            end
        end
        % Se il rapporto è già stato calcolato in un passo precedente, lo manteniamo invariato
    end
end
