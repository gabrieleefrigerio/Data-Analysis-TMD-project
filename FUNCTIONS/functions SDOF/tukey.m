function [HAMM, FIRST] = tukey(hammer, L, flag)
% INPUTS:
% 1) hammer: matrix containing the time histories of the hammer
% hits
% 2) L: width of the window espressed as number of samples
% (suggested: fsamp/100)
% 3) Flag: 1 if we want to window the signal, another number otherwise

% OUTPUTS
% 1) HAMM: matrix containting the windowed signals of the hammer
% 2) FIRST: index indicating the first sample of all the windows
% (this will be used for cutting the signal)

N = size(hammer,1);
wind = zeros(N,size(hammer,2));

if flag == 1
    FIRST = 1e10;   % initialized with a high value

    for ii = 1:size(hammer,2)
        idx_max_hamm = find(max(hammer(:,ii))==hammer(:,ii));
        % generates a Tukey window properly centered (0.8 con be kept fixed)
        wind(idx_max_hamm-floor((L-1)/2):idx_max_hamm+ceil((L-1)/2),ii) = tukeywin(L,0.2);
        HAMM(:,ii) = wind(:,ii).*hammer(:,ii);  % windowed signal
        if idx_max_hamm-floor((L-1)/2) < FIRST
            FIRST = idx_max_hamm-floor((L-1)/2);
        end
    end
else
    HAMM = hammer;
    FIRST = 1;
end

% plots to verify the windowning
for ii = 1:size(hammer,2)
    idx_max_hamm = find(max(hammer(:,ii))==hammer(:,ii));

    % figure
    % plot(hammer(:,ii)); hold on; grid on
    % plot(wind(:,ii)*max(hammer(:,ii)))
    % title(['Test ', num2str(ii)])
    % xlim([idx_max_hamm-(L-1)/2,idx_max_hamm+(L-1)/2])
end

end
