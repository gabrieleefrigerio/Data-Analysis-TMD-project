function [FRF, COH, SXX, SYY] = get_frf(data,hammer)
% Computes the frequency response functions using the esimators, computes
% the coherence
% INPUTS:
% 1) data: 3D matrix from go_to_mat
% 2) hammer: 2D matrix from tukey

% OUTPUTS:
% 1) FRF: matrix N/2 x M containing the frequency response function of
% each accelerometer
% 2) COH: matrix N/2 x M containing the coherence of each accelerometer
% 3) SXX: matrix N/2 x P containing the autospectrum of the hammer signal
% of each test
% 4) SYY: matrix N/2 x M x P containing the autospectrum of each
% accelerometer of each test
% N: number of samples
% M: number of accelerometers
% P: number of tests


N = size(data,1);

%preallocating vectors
Sxx = zeros(N, size(hammer, 2));
Syy = zeros(N, size(data, 2), size(data, 3));
Sxy = zeros(N, size(data, 2), size(data, 3));
Syx = zeros(N, size(data, 2), size(data, 3));

for ii = 1:size(data,2)
    for jj = 1:size(data,3)
        % fast Fourier transform
        x = fft(hammer(:,jj))/N;
        y = fft(data(:,ii,jj))/N;

        % auto and cross spectra (double sided)
        Sxx(:,jj)    = 4*conj(x).*x;
        Syy(:,ii,jj) = 4*conj(y).*y;
        Sxy(:,ii,jj) = 4*conj(x).*y;
        Syx(:,ii,jj) = 4*conj(y).*x;
    end
end

% single sided spectra
Sxx = Sxx(1:floor(N/2),:);
Sxy = Sxy(1:floor(N/2),:,:);
Syy = Syy(1:floor(N/2),:,:);
Syx = Syx(1:floor(N/2),:,:);

% H1 and H2 estimators
H1 = mean(Sxy,3)./mean(Sxx,2);  
H2 = mean(Syy,3)./mean(Syx,3); 

COH = abs(H1./H2);  % coherence
FRF = H1;           % choice of the frf between H1 and H2
SXX = Sxx;
SYY = Syy;

end