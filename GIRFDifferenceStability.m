%% 3dB bandwidth analysis

%% File name and path
resultPath = "../ISMRM2022Results/";
resultFileName2 = "2021Jun_Gx.mat";
resultFileName1 = "2020Nov_Gx.mat";
% resultFileName2 = "2021Jun_Gy.mat";
% resultFileName1 = "2020Nov_Gy.mat";
% resultFileName2 = "2021Jun_Gz.mat";
% resultFileName1 = "2020Nov_Gz.mat";

% This will load the following variables:
% GIRF_FT, b0ec_FT (if have), dwellTimeSig, isAvgRepetition, roPts, roTime
load(strcat(resultPath, resultFileName1));
GIRF_FT1 = GIRF_FT; clear GIRF_FT;
load(strcat(resultPath, resultFileName2));
GIRF_FT2 = GIRF_FT; clear GIRF_FT;

%% Calculation
threashold = 0.8; % Threshold of cut-off rate

freq_fullrange = 1 / (dwellTimeSig / 1e6) / 1e3; % Full spectrum width, in unit of kHz
freq = linspace(-freq_fullrange/2, freq_fullrange/2, roPts);
freq = freq(:);

dispFreqRange = [-30, 30]; % in unit of kHz

% [~, freqDispIndexStart] = min(abs(freq - dispFreqRange(1)));
% [~, freqDispIndexEnd] = min(abs(freq - dispFreqRange(end)));
% freqDispIndex = freqDispIndexStart : freqDispIndexEnd;
% freqDispIndex = freqDispIndex(:);

GIRF_FT_mean1 = abs(mean(GIRF_FT1,2));
GIRF_FT_mean2 = abs(mean(GIRF_FT2,2));
% GIRF_FT_mean1 = GIRF_FT_mean1(freqDispIndex);
% GIRF_FT_mean2 = GIRF_FT_mean2(freqDispIndex);

% % Max amplitude and indexes for amplitude with 3dB decay
% maxGIRF1 = max(GIRF_FT_mean1); % Need to use center points or the high frequency part with lowest SNR will affect this level
% maxGIRF2 = max(GIRF_FT_mean2);
% [~,indexLeft1] = min(abs(GIRF_FT_mean1(1:end/2) - maxGIRF1*threashold));
% [~,indexRight1] = min(abs(GIRF_FT_mean1(end/2+1 : end) - maxGIRF1*threashold));
% indexRight1 = indexRight1 + length(GIRF_FT_mean1)/2;
% [~,indexLeft2] = min(abs(GIRF_FT_mean2(1:end/2) - maxGIRF2*threashold));
% [~,indexRight2] = min(abs(GIRF_FT_mean2(end/2+1 : end) - maxGIRF2*threashold));
% indexRight2 = indexRight2 + length(GIRF_FT_mean2)/2;

% freqDist1 = freq(freqDispIndex(indexRight1)) - freq(freqDispIndex(indexLeft1))
% freqDist2 = freq(freqDispIndex(indexRight2)) - freq(freqDispIndex(indexLeft2))

freq1 = 3.2; % Freq to observe according to EPI/Spiral Traj, in kHz
freq2 = 14.6;
[~,indexLeft1] = min(abs(freq + freq1));
[~,indexRight1] = min(abs(freq - freq1));
[~,indexLeft2] = min(abs(freq + freq2));
[~,indexRight2] = min(abs(freq - freq2));

%% Plot the difference of GIRF
figure(222);
clf;
set(gcf,'color','white');
plot(freq, abs(GIRF_FT_mean1), 'r', 'LineWidth', 2);
hold on;
plot(freq, abs(GIRF_FT_mean2), 'k', 'LineWidth', 2);
hold on;
plot(freq, abs(GIRF_FT_mean1 - GIRF_FT_mean2), 'b', 'LineWidth', 1);
hold on;
xline(freq(indexLeft1), '--', 'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880]);
hold on;
xline(freq(indexRight1), '--', 'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880]);
hold on;
xline(freq(indexLeft2), '--', 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250]);
hold on;
xline(freq(indexRight2), '--', 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250]);
% plot([freq(freqDispIndex(indexLeft1)), freq(freqDispIndex(indexLeft1))], [0, 1.1], 'Color', 'r');
% hold on;
% plot([freq(freqDispIndex(indexRight1)), freq(freqDispIndex(indexRight1))], [0, 1.1], 'Color', 'r');
% hold on;
% plot([freq(freqDispIndex(indexLeft2)), freq(freqDispIndex(indexLeft2))], [0, 1.1], 'Color', 'k');
% hold on;
% plot([freq(freqDispIndex(indexRight2)), freq(freqDispIndex(indexRight2))], [0, 1.1], 'Color', 'k');
xlim(dispFreqRange); ylim([0, 1.1]);
xlabel('Frequency [kHz]','FontSize', 18); ylabel('Magnitude of GIRF','FontSize', 18);
% title('GIRF in Frequency Domain', 'FontSize', 26);
% legend('Initial', '8 Months Later', 'Difference','FontSize', 18);

girfDiffLeftFreq1 = GIRF_FT_mean1(indexLeft1) - GIRF_FT_mean2(indexLeft1)
girfDiffRightFreq1 = GIRF_FT_mean1(indexRight1) - GIRF_FT_mean2(indexRight1)
girfDiffLeftFreq2 = GIRF_FT_mean1(indexLeft2) - GIRF_FT_mean2(indexLeft2)
girfDiffRightFreq2 = GIRF_FT_mean1(indexRight2) - GIRF_FT_mean2(indexRight2)

girfRelDiff = abs(GIRF_FT_mean1 - GIRF_FT_mean2) ./ GIRF_FT_mean1;
girfRelDiffFreq1 = max(girfRelDiff(indexLeft1:indexRight1))
girfRelDiffFreq2 = max(girfRelDiff(indexLeft2:indexRight2))

%%
% %% Plot B0 Eddy Current (k0) in Frequency Domain
% freq_fullrange = 1 / (dwellTimeSig / 1e6) / 1e3; % Full spectrum width, in unit of kHz
% freq = linspace(-freq_fullrange/2, freq_fullrange/2, roPts);
% 
% dispFreqRange = [-30, 30]; % in unit of kHz
% 
% % freq_index = 4551:5451; % -9kHz ~ +9kHz
% % freq_index = 4501:5501; % -10kHz ~ +10kHz
% % freq_index = 4251:5751; % -15kHz ~ +15kHz
% % freq_index = 4001:6001; % -20kHz ~ +20kHz
% % freq_index = 3501:6501; % -30kHz ~ +30kHz
% % freq_index = 3001:7001; % -40kHz ~ +40kHz
% 
% figure(888);
% set(gcf,'color','white');
% 
% % subplot(211);
% if isAvgRepetition == 0
%     plot(freq, abs(mean(b0ec_FT,2)));
% else
%     plot(freq, abs(b0ec_FT));
% end
% xlim(dispFreqRange); ylim([0, 500]);
% xlabel('Frequency [kHz]', 'FontSize', 14); ylabel('Magnitude of B0 EC [AU]', 'FontSize', 14);
% % title('Magnitude of B0 Eddy Current in Frequency Domain', 'FontSize', 18);
