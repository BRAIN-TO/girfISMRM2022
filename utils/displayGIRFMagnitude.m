function displayGIRFMagnitude(GIRF_FT, fullFreqRange, dispFreqRange, figNum)
%Function for displaying GIRF calculation results
% function displayGIRFMagnitude(GIRF_FT, fullFreqRange, dispFreqRange, figNum)
% Input arguments:
%   GIRF_FT: GIRF in frequency domain, column vector.
%   fullFreqRange: points covering full frequency range, column vector.
%   dispFreqRange: min and max values of displaying frequency range, 1 × 2 vector. (optional)
%   figNum: figure number. (optional)
%    
% Author: Zhe "Tim" Wu
% Created: March 10, 2022

    if nargin == 2
        dispFreqRange = [min(fullFreqRange), max(fullFreqRange)];
        figure;
    elseif nargin == 3
        figure;
    elseif nargin == 4
        figure(figNum);
    else
        error('Invalid input arguments, check help document of the function.');
    end
        
    set(gcf,'color','white');

    GIRF_FT_mean = mean(GIRF_FT,2);
    GIRF_FT_std = std(GIRF_FT, 0, 2);
    GIRF_FT_mean_abs = abs(GIRF_FT_mean);
    GIRF_FT_mean_phase = angle(GIRF_FT_mean);
    GIRF_FT_std_abs = abs(GIRF_FT_std);
    
    % Display magnitude (mean and standard deviation)
    subplot(211);
    errorBarColor = [255 99 71]/255; % Error bar in Tomato Red Color
    patch = fill([fullFreqRange',fliplr(fullFreqRange')], [GIRF_FT_mean_abs' - GIRF_FT_std_abs', fliplr(GIRF_FT_mean_abs' + GIRF_FT_std_abs')], errorBarColor);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0.8);
    hold on;  
    plot(fullFreqRange, GIRF_FT_mean_abs, 'b', 'LineWidth', 1);
    xlim(dispFreqRange); ylim([0, 1.1]);
    xlabel('Frequency [kHz]', 'FontSize', 14); ylabel('Magnitude of GIRF [AU]', 'FontSize', 14);
    title('Magnitude of GIRF in Frequency Domain', 'FontSize', 14);
    hold off;
    
    % Fitting phase part for gradient temporal delay
    fitRange = [-9, 9]; % in kHz
    [~, freqFitIndexStart] = min(abs(fullFreqRange - fitRange(1)));
    [~, freqFitIndexEnd] = min(abs(fullFreqRange - fitRange(end)));
    freqFitIndex = freqFitIndexStart : freqFitIndexEnd;
    freqFitIndex = freqFitIndex(:);
    fitX = fullFreqRange(freqFitIndex);
    fitY = GIRF_FT_mean_phase(freqFitIndex);
    
    p1 = polyfit(fitX, fitY, 1);
    temporalDelay1 = p1(1)/1e3/(2*pi) * 1e6; % in unit of us
    disp(['Delay:', num2str(temporalDelay1), ' us']);

    % Display phase (mean only)
    subplot(212);
    plot(fullFreqRange, GIRF_FT_mean_phase, 'b', 'LineWidth', 1);
    xlim(dispFreqRange); ylim([-4,4]);
    xlabel('Frequency [kHz]', 'FontSize', 14); ylabel('Phase of GIRF [rad]', 'FontSize', 14);
    title(strcat('Phase of GIRF in Frequency Domain, delay is ', num2str(temporalDelay1), ' us'), 'FontSize', 14);

    

end