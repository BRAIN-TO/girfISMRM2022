%% Load the files
gradPath = '../ISMRM2022Results/';
% gradFileName = 'FigTrajGradSlew_311.mat'; % EPI
% titleStr = 'Frequency Spectrum For EPI Readout'; % EPI
gradFileName = 'FigTrajGradSlew_511.mat'; % Spiral
titleStr = 'Frequency Spectrum For Spiral Readout'; % Spiral

load(strcat(gradPath, gradFileName)); % will load struct "data"

% struct "data" contains following domains:
% FT_G, freq_kHz, timeGradient, timeSlew, time, k, g, s

% figure; plot(data.timeGradient, data.g);
% figure; plot(data.k(:,1), data.k(:,2));
% figure; plot(data.freq_kHz, abs(data.FT_G));

startIndex = int16(length(data.freq_kHz) / 2);
freqDisp = data.freq_kHz(startIndex:end);
integralAmp = cumsum(abs(data.FT_G(startIndex:end,:)).^2,1);
integralAmp(:,1) = integralAmp(:,1) / integralAmp(end,1);
integralAmp(:,2) = integralAmp(:,2) / integralAmp(end,2);
[~, xIndex99] = min(abs(integralAmp(:,1) - 0.99)); % Find threshold for 99% power accumulation
[~, yIndex99] = min(abs(integralAmp(:,2) - 0.99));
disp(['99% Power Threshold on Gx: ', num2str(freqDisp(xIndex99)), 'kHz']); 
disp(['99% Power Threshold on Gy: ', num2str(freqDisp(yIndex99)), 'kHz']);

% figure; plot(data.freq_kHz(startIndex:end), integralAmp, 'LineWidth', 2);
figure;
set(gcf,'color','white');
plot(freqDisp, abs(data.FT_G(startIndex:end,1)), 'LineWidth', 3, 'Color', 'r');
hold on;
plot(freqDisp, abs(data.FT_G(startIndex:end,2)), 'LineWidth', 3, 'Color', 'k');
hold on;
xline(freqDisp(xIndex99), '--r', 'LineWidth', 2);
hold on;
xline(freqDisp(yIndex99), '--k', 'LineWidth', 2);
title(titleStr, 'FontSize', 18);
xlabel('Frequency [kHz]', 'FontSize', 14); ylabel('Amplitude [AU]', 'FontSize', 14);
legend('|FT(Gx)|', '|FT(Gy)|', '99% Gx Power', '99% Gy Power', 'FontSize', 14);

xlim([0, 5]); % For Spiral 511, in unit of kHz
% xlim([0, 20]); ylim([0, 7000]); % For EPI 311 , in unit of kHz