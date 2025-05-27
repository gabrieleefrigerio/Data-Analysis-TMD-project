function [DATA, HAMMER, TIME, FVEC] = cut_data(data, hammer, time, first_samp)
% Cuts the data according to the index from the function tukey 
% INPUTS:
% 1) data: 3D matrix from go_to_mat
% 2) hammer: 2D matrix from tukey
% 3) time: time vector
% 4) first_samp: index from tukey
% 5) fsamp: sampling frequency
%
% OUTPUTS:
% 1) DATA: cut 3D matrix
% 2) HAMMER: cut 2D matrix
% 3) TIME: cut time vector
% 4) FVEC: frequency vector

DATA = data(first_samp:end,:,:);
TIME = time(first_samp:end,:) - time(first_samp);
HAMMER = hammer(first_samp:end,:);
df = 1/TIME(end);
N = size(DATA,1);
FVEC = transpose(0:df:(floor(N/2)-1)*df);

% plots to verify the cut
% for ii = 1:size(DATA,3)
%     figure
%     plot(TIME, DATA(:,:,ii)); hold on; grid on
%     title(['Test ', num2str(ii)])
% end

end