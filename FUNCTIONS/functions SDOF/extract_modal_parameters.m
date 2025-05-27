function [freq, damping] = extract_modal_parameters(A, dt)
    % Autovalori
    lambda = eig(A);
    
    % Modal frequencies and damping from eigenvalues
    s = log(lambda) / dt;  % continuo
    freq = abs(imag(s)) / (2*pi);  % Hz
    damping = -real(s) ./ abs(s);  % rapporto smorzamento
end
