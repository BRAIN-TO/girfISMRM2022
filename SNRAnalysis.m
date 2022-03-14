%% SNR analysis

%% File name and path
resultPath = "../ISMRM2022Results/";
% resultFileName1 = "2021Jun_Gx.mat";
% resultFileName2 = "2021Jun_PosNeg_Gx.mat";
% resultFileName1 = "2021Jun_Gy.mat";
% resultFileName2 = "2021Jun_PosNeg_Gy.mat";
resultFileName1 = "2021Jun_Gz.mat";
resultFileName2 = "2021Jun_PosNeg_Gz.mat";

% This will load the following variables:
% GIRF_FT, b0ec_FT (if have), dwellTimeSig, isAvgRepetition, roPts, roTime
load(strcat(resultPath, resultFileName1));
GIRF_FT1 = GIRF_FT; clear GIRF_FT;
load(strcat(resultPath, resultFileName2));
GIRF_FT2 = GIRF_FT; clear GIRF_FT;

%% SNR Calculation
GIRF_FT_mean1 = abs(mean(GIRF_FT1,2));
GIRF_FT_std1 = abs(std(GIRF_FT1, 0, 2));
SNR1 = GIRF_FT_mean1 ./ GIRF_FT_std1;

GIRF_FT_mean2 = abs(mean(GIRF_FT2,2));
GIRF_FT_std2 = abs(std(GIRF_FT2, 0, 2));
SNR2 = GIRF_FT_mean2 ./ GIRF_FT_std2;

%% SNR sliding window smoothing
windowLength = 100; % in pts
SNRSmooth1 = movmean(SNR1, windowLength);
SNRSmooth2 = movmean(SNR2, windowLength);

maxSNRRatio = max(SNRSmooth2) ./ max(SNRSmooth1)

%% Plot SNR
freq_fullrange = 1 / (dwellTimeSig / 1e6) / 1e3; % Full spectrum width, in unit of kHz
freq = linspace(-freq_fullrange/2, freq_fullrange/2, roPts);
freq = freq(:);

dispFreqRange = [-30, 30]; % in unit of kHz

figure(111);
set(gcf,'color','white');
% plot(freq, SNR1, 'b', 'LineWidth', 1);
plot(freq, SNRSmooth1, 'r', 'LineWidth', 1);
hold on;
% plot(freq, SNR2, 'r', 'LineWidth', 1);
plot(freq, SNRSmooth2, 'k', 'LineWidth', 1);
hold on;
xlim(dispFreqRange);
xlabel('Frequency [kHz]','FontSize', 14); ylabel('SNR [AU]','FontSize', 14);
title('SNR of GIRF in Frequency Domain','FontSize', 18);
hold off;
legend('Origin', 'Optimized', 'FontSize', 14);

figure(112);
set(gcf,'color','white');
plot(freq, SNRSmooth2./SNRSmooth1, 'b', 'LineWidth', 1);
xlim(dispFreqRange);
xlabel('Frequency [kHz]','FontSize', 14); ylabel('SNR [AU]','FontSize', 14);
title('SNR of GIRF in Frequency Domain','FontSize', 18);
hold off;

% figure(113);
% powerbw(fftshift(fft(fftshift(SNRSmooth1,1),[],1),1), 1e6/5/(10/3));

% Mean of SNR ratio between display range
% Noisy part < -30 kHz or > 30 kHz was removed
[~, freqDispIndexStart] = min(abs(freq - dispFreqRange(1)));
[~, freqDispIndexEnd] = min(abs(freq - dispFreqRange(end)));
freqDispIndex = freqDispIndexStart : freqDispIndexEnd;
freqDispIndex = freqDispIndex(:);
meanSNRRatio = mean(SNRSmooth2(freqDispIndex) ./ SNRSmooth1(freqDispIndex))

