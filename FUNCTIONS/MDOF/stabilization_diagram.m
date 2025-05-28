function [modal_results] = stabilization_diagram(IRF, K, max_order, FRF)
% stabilization_diagram - Build stabilization diagram and extract modal parameters
%
% Inputs:
%   IRF       - Struct containing impulse response function data
%   K         - Parameter m+n defining Hankel matrix size
%   max_order - Maximum model order to test
%   FRF       - Struct with frequency response function data
%
% Output:
%   modal_results - Struct with modal frequencies, modes, and poles

    % Initialize cell arrays to store frequencies, damping, orders, and modes
    all_freqs = {};
    all_zeta = {};
    all_orders = {};
    all_modes = {};  % Full eigenvectors (mode shapes)

    figure; hold on;

    % Plot FRF on left y-axis (log scale)
    yyaxis left
    semilogy(FRF.freq, FRF.FRF, 'b-', 'LineWidth', 2);
    ylabel('FRF magnitude');

    % Extract IRF data dimensions
    h = IRF.irf;  % Nt x nChan x nAcq
    [Nt, nChan, nAcq] = size(h);

    % Switch to right y-axis for modal data plotting
    yyaxis right
    dt = IRF.time(2) - IRF.time(1);

    % Loop over model orders from 1 to max_order
    for i = 1:(max_order)
        m = 100+i;
        n = K - m;

        % Build Hankel matrices (try-catch to skip invalid orders)
        try
            [H0, H1] = build_hankel(IRF, m, n);
        catch
            warning("⚠️  Skipping m=%d, n=%d: IRF too short\n", m, n);
            continue;
        end

        %size_eig = 2 * min(m, n);
        size_eig = 2*i;

        % Check rank condition on H0 matrix
        if rank(H0) < size_eig
            warning("⚠️  Order too high for rank: H0 rank = %d\n", rank(H0));
            continue;
        end

        % Compute truncated SVD of H0
        [U, ~, ~] = svds(H0, size_eig);
        Ur = U;

        % Compute state matrix A
        H0_prime = pinv(Ur' * H0);
        A = (Ur' * H1) * H0_prime;

        % Eigen decomposition of A
        [V, D] = eig(A);  % V = reduced eigenvectors, D = eigenvalues

        % Calculate continuous-time poles s = log(lambda)/dt
        s = log(diag(D)) / dt;

        % Extract modal frequencies (Hz) and damping ratios
        freq_Hz = abs(imag(s)) / (2 * pi);
        zeta = -real(s) ./ abs(s);

        % Calculate full mode shapes by expanding reduced eigenvectors
        Phi_wave = Ur * V;
        Phi = Phi_wave(1:nChan, :);

        % Store results in cell arrays indexed by order (size_eig/2 = model order)
        all_freqs{size_eig / 2} = freq_Hz(:);
        all_zeta{size_eig / 2} = zeta(:);
        all_orders{size_eig / 2} = repmat(size_eig / 2, length(freq_Hz), 1);
        all_modes{size_eig / 2} = num2cell(Phi, 1); % store each mode shape as a separate cell
    end

    %% Stability Analysis by Frequency Binning

    % Use last order's frequencies as bin centers
    ref_freqs = all_freqs{end};
    bin_centers = ref_freqs(:);

    % Define bin widths as 0.1% of each frequency
    bin_widths = 0.001 * bin_centers;
    n_bins = length(bin_centers);

    % Initialize bin counters
    bin_counts = zeros(n_bins, 1);

    % Count frequencies falling into each bin across all orders
    for i = 1:length(all_freqs)
        freqs = all_freqs{i};
        if isempty(freqs)
            continue;
        end
        for f = reshape(freqs, 1, [])
            for b = 1:n_bins
                half_width = bin_widths(b) / 2;
                if f >= (bin_centers(b) - half_width) && f <= (bin_centers(b) + half_width)
                    bin_counts(b) = bin_counts(b) + 1;
                    break; % frequency belongs to only one bin
                end
            end
        end
    end

    % Identify stable frequency bins (appearing >= 10 times)
    stable_bins = bin_counts >= 10;
    stable_freqs = bin_centers(stable_bins);

    % Plot stable frequencies on right axis
    yyaxis right
    plot(bin_centers(stable_bins), zeros(sum(stable_bins), 1), 'go', 'MarkerSize', 8, ...
        'DisplayName', 'Stable modes');

    legend show;

    %% Stability Analysis on Damping (zeta)

    ref_zeta = all_zeta{end};
    zeta_centers = ref_zeta(:);
    zeta_widths = 0.05 * zeta_centers; % 5% damping tolerance
    n_zbins = length(zeta_centers);
    zeta_counts = zeros(n_zbins, 1);

    % Count damping ratios in each bin
    for i = 1:length(all_zeta)
        zetas = all_zeta{i};
        if isempty(zetas)
            continue;
        end
        for z = reshape(zetas, 1, [])
            for b = 1:n_zbins
                half_width = zeta_widths(b) / 2;
                if z >= (zeta_centers(b) - half_width) && z <= (zeta_centers(b) + half_width)
                    zeta_counts(b) = zeta_counts(b) + 1;
                    break;
                end
            end
        end
    end

    % Identify stable damping bins
    stable_zeta_bins = zeta_counts >= 10;
    stable_zeta_values = zeta_centers(stable_zeta_bins);

    % Filter stable frequencies by stable damping
    stable_freqs_filtered = [];
    for i = 1:length(stable_freqs)
        f = stable_freqs(i);

        % Find nearest bin index
        [~, idx] = min(abs(bin_centers - f));
        corresponding_zeta = ref_zeta(idx);

        % Check if damping is within stable damping bins
        is_stable_zeta = false;
        for j = 1:length(stable_zeta_values)
            center = stable_zeta_values(j);
            width = 0.05 * center / 2;
            if corresponding_zeta >= (center - width) && corresponding_zeta <= (center + width)
                is_stable_zeta = true;
                break;
            end
        end

        if is_stable_zeta
            stable_freqs_filtered(end+1, 1) = f; %#ok<AGROW>
        end
    end

    %% Modal Assurance Criterion (MAC) Filtering

    stable_freqs_mac_filtered = [];

    for i = 1:length(stable_freqs_filtered)
        f = stable_freqs_filtered(i);

        % Find nearest frequency bin
        [~, idx] = min(abs(bin_centers - f));

        % Reference mode shape for this frequency
        phi_ref = all_modes{end}{idx};
        if isempty(phi_ref)
            continue;
        end

        macs = [];

        % Compare mode shape with neighbors within ±0.05% frequency
        freq_bin_width = 0.001 * f / 2;
        for j = 1:length(bin_centers)
            fj = bin_centers(j);
            if abs(fj - f) <= freq_bin_width
                phi = all_modes{end}{j};
                if isempty(phi)
                    continue;
                end

                % Calculate MAC value
                mac = abs(phi_ref' * phi)^2 / ((phi_ref' * phi_ref) * (phi' * phi));
                macs(end+1) = mac;
            end
        end

        % Keep frequency if MAC ≥ 0.8 with any neighbor
        if any(macs >= 0.8)
            stable_freqs_mac_filtered(end+1, 1) = f; %#ok<AGROW>
        end
    end

    %% Extract selected modes and poles corresponding to stable frequencies

    ref_modes = all_modes{end};
    ref_poles = s;  % poles from last iteration

    selected_modes = {};
    selected_poles = [];

    for i = 1:length(stable_freqs_mac_filtered)
        f = stable_freqs_mac_filtered(i);

        [~, idx] = min(abs(ref_freqs - f));
        phi = ref_modes{idx};
        if ~isempty(phi)
            selected_modes{end+1, 1} = phi;
            selected_poles(end+1, 1) = ref_poles(idx);
        end
    end

    % Package modal results in output struct
    modal_results = struct;
    modal_results.eigenfreq = stable_freqs_mac_filtered;
    
    nModes = length(selected_modes);          % Number of modes
    nOut = length(selected_modes{1});         % Number of channels (size of mode vector)
    
    modal_results.modes = zeros(nOut, nModes);  % Preallocate modes matrix
    for k = 1:nModes
        modal_results.modes(:, k) = selected_modes{k};
    end
    
    modal_results.poles = selected_poles;
    modal_results.damping = -real(selected_poles) ./ abs(selected_poles);

    %% Plot stabilization diagram details

    hold on; grid on;
    title('Stabilization Diagram');
    xlabel('Frequency [Hz]');
    ylabel('Model order');

    % Plot FRF magnitude on right y-axis
    yyaxis right
    plot(FRF.freq, abs(FRF.FRF), 'b-', 'LineWidth', 1.5);
    ylabel('FRF Magnitude');

    % Plot estimated modes for each order with transparency
    colors = lines(max_order);
    for order = 1:max_order
        idx = order;
        if idx <= length(all_freqs) && ~isempty(all_freqs{idx})
            scatter(all_freqs{idx}, order * ones(size(all_freqs{idx})), 20, ...
                colors(order,:), 'filled', 'MarkerFaceAlpha', 0.3);
        end
    end

    % Highlight final stable modes with red pentagram markers
    stable_orders = [];
    for i = 1:length(stable_freqs_mac_filtered)
        f = stable_freqs_mac_filtered(i);
        ord_idx = find(cellfun(@(c) any(abs(c - f) < 1e-6), all_freqs), 1);
        if isempty(ord_idx)
            ord_idx = max_order;
        end
        stable_orders(end+1) = ord_idx; %#ok<AGROW>
        scatter(f, ord_idx, 80, 'r', 'p', 'filled', 'MarkerEdgeColor', 'k');
    end

    legend({'FRF Magnitude', 'Modi per ordine', 'Modi stabili finali'}, 'Location', 'best');
hold off;

end