function H0_pinv = compute_pseudoinverse(H0, tol)
% compute_pseudoinverse Calcola la pseudoinversa di H0 via SVD
%
%   H0_pinv = compute_pseudoinverse(H0, tol)
%   tol: soglia per il troncamento dei valori singolari (opzionale)

    if nargin < 2
        tol = max(size(H0)) * eps(norm(H0));  % default numerico
    end

    [U, S, V] = svd(H0, 'econ');

    % Inverti solo i valori singolari significativi
    s = diag(S);
    r = sum(s > tol);  % rango numerico

    if r == 0
        warning('Pseudoinversa mal condizionata: rango nullo.');
        H0_pinv = zeros(size(H0'));
        return;
    end

    S_inv = diag(1 ./ s(1:r));
    H0_pinv = V(:,1:r) * S_inv * U(:,1:r)';  % H0^†
end
