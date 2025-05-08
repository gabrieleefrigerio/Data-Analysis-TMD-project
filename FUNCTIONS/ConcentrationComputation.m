function [data] = ConcentrationComputation(data)
    % ConcentrationComputation - Calcola il rapporto di concentrazione di martensite 
    % in un materiale a memoria di forma (SMA) in base alla temperatura attuale e 
    % a quella del passo precedente. Il rapporto varia tra 1 (solo martensite) e 
    % 0 (solo austenite), con interpolazione lineare nelle regioni di transizione.
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
        if data.T < data.As
            data.ratio = 1;
        elseif data.T < data.Af
            data.ratio = (data.Af - data.T) / (data.Af - data.As);
        else
            data.ratio = 0;
        end

    elseif data.T < data.T_prev
        % --- Raffreddamento ---
        if data.T < data.Mf
            data.ratio = 1;
        elseif data.T < data.Ms
            data.ratio = (data.Mf - data.T) / (data.Mf - data.Ms);
        else
            data.ratio = 0;
        end

    else
        % --- Temperatura costante: corpo a regime termico ---
        if ~isfield(data, 'ratio')
            if data.T <= data.Mf
                data.ratio = 1;
            elseif data.T >= data.Af
                data.ratio = 0;
            elseif data.T > data.Mf && data.T < data.Ms
                % Zona di transizione martensitica → interpolazione
                data.ratio = (data.Mf - data.T) / (data.Mf - data.Ms);
            elseif data.T > data.As && data.T < data.Af
                % Zona di transizione austenitica → interpolazione
                data.ratio = (data.Af - data.T) / (data.Af - data.As);
            else
                % Tra Ms e As, regione non specificata chiaramente
                data.ratio = 0.5;
            end
        end
    end
end
