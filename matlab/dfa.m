function [scales, fluctuation, alpha] = dfa(x, scales, order, showPlot)
    % Perform Detrended Fluctuation Analysis on x
    % Parameters:
    %   x (column vector): 1D numeric array containing time series x
    %   scales (numeric array): Array of scales to calculate fluctuations
    %   order (integer): Order of polynomial fit (1 = linear fit)
    %   showPlot (logical): draw the log-log fit. Default false, so that
    %       batch and cluster runs do not create a figure they cannot close.
    %
    % Returns:
    %   scales: The scales that were entered as input
    %   fluctuations: Variability measured at each scale with RMS
    %   alpha value: Value quantifying the variability in the time series
    %
    % References:
    %   Peng, C. K., Havlin, S., Stanley, H. E., & Goldberger, A. L. (1995).
    %   Quantification of scaling exponents and crossover phenomena in 
    %   nonstationary heartbeat time series. Chaos: An Interdisciplinary 
    %   Journal of Nonlinear Science, 5(1), 82-87.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

%% ========================================================================
%                          ------ EXAMPLE ------

%      - Generate random x
%      x = randn(5000,1); 
      
%      - Create a vector of the scales you want to use
%      scales = [10, 20, 40, 80, 160, 320, 640, 1280, 2560];
      
%      - Set a detrending order. Use 1 for a linear detrend.
%      order = 1;
      
%      - run dfa function
%      [s, f, a] = dfa(x, scales, order, 1)

%% ========================================================================

% Check if the x is a column vector and if not, transpose to make it
% one.
if nargin < 4 || isempty(showPlot)
    showPlot = false;
end

if size(x, 2) > size(x, 1) % x should be column vector
    x = x';
end


    % Integrate the x
    integrated_data = cumsum(x - mean(x));

    fluctuation = zeros(size(scales));

    for idx = 1:length(scales)
        scale = scales(idx);

        % Divide x into non-overlapping chunks of size 'scale'
        chunks = floor(length(x) / scale);
        %disp(chunks)
        ms = 0;

        for i = 1:chunks
            chunk_start = (i - 1) * scale + 1;
            chunk_end = i * scale;
            this_chunk = integrated_data(chunk_start:chunk_end);
            chunk_idx = 1:length(this_chunk);

            % Fit polynomial (default is linear, i.e., order=1)
            coeffs = polyfit(chunk_idx, this_chunk, order);
            fit = polyval(coeffs, chunk_idx);

            % Detrend and calculate RMS for this chunk
            ms = ms + mean((this_chunk' - fit).^2);
        end

        % Calculate average RMS for this scale
        fluctuation(idx) = sqrt(ms./chunks);
    end

    % Perform linear regression on the log-log x
    log_scales = log(scales);
    log_fluctuation = log(fluctuation);

    % Check for NaN values in log_fluctuation and remove them
    nan_indices = isnan(log_fluctuation);
    log_scales(nan_indices) = [];
    log_fluctuation(nan_indices) = [];

    p = polyfit(log_scales, log_fluctuation, 1);
    alpha = p(1);

    % Calculate R-squared value
    y_fit = polyval(p, log_scales);
    ssr = sum((log_fluctuation - y_fit).^2);
    sst = sum((log_fluctuation - mean(log_fluctuation)).^2);
    rsquared = 1 - ssr / sst;


    if showPlot
        % Plot scales vs. fluctuation values
        loglog(scales, fluctuation, 'o-k', MarkerFaceColor='red', MarkerSize = 8, LineWidth= 1.5);
        hold on;
        xlabel('Scale (log)');
        ylabel('Fluctuation (log)');

        % Add alpha and R-squared as text labels
        str = ['Alpha = ', num2str(alpha, '%.3f') newline 'R^2 = ', num2str(rsquared, '%.3f')];
        dim = [.7 .05 .2 .2]; % textbox location
        annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineWidth', 1, 'BackgroundColor', 'white');


        % Add the regression line to the plot
        loglog(scales, exp(polyval(p, log_scales)), 'Color', '#4DBEEE', 'LineStyle', '--', 'LineWidth', 1.5);

        title('Detrended Fluctuation Analysis');
        grid on;
        grid minor;
    end
end
