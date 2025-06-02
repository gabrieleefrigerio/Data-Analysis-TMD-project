function [modal_mass] = find_modal_mass(FRF,modal_results)
    % Dati di input
    Hpq = FRF.FRF;                         % (Nf x 1)
    psi = modal_results.modes;             % (n x 2N)
    sr = modal_results.poles;              % (2N x 1)
    omega = FRF.freq *2*pi;                    % (Nf x 1)
    xi = modal_results.damping;
    wr = modal_results.eigenfreq*2*pi;
    % Assumiamo p e q noti, ad esempio:
    p = 1;
    q = 1;
    
    % Estrai i vettori colonna per psi_p e psi_q
    psi_p = psi(p, :).';                 % (2N x 1)
    psi_q = psi(q, :).';                 % (2N x 1)
    
    % Calcolo della matrice A
    Nf = length(omega);
    num_poli = length(sr);
    
    A = zeros(Nf, num_poli + 2);         % ultima colonna per 1/omega^2, penultima per 1
    
    for k = 2:Nf
        jw = 1j * omega(k);
        for i = 1:num_poli
            A(k, i) = (psi_p(i) * psi_q(i)) / (jw - sr(i));
        end
        A(k, end-1) = 1;
        A(k, end) = 1 / omega(k)^2;
    end
    
    % Ora Hpq ≈ A * x, con x = [Q_1, ..., Q_2N, R_U, R_L].'
    % Puoi risolverlo per x con:
    x = A \ Hpq;
    
    % x conterrà: [Q; RU; RL]
    modal_mass = zeros(num_poli, 1);
    for ii = 1:num_poli
        Qr = x(ii);
        modal_mass(ii) = 1 ./ (2j .* Qr.* wr(ii) .* sqrt(1 - xi(ii).^2));
        modal_mass(ii) = abs(modal_mass(ii));
    end

end