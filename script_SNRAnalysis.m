% script_SNRAnalysis.m
%Script for GIRF SNR analysis, generating subfigures on the second row of Figure 4
% To analyze the SNR difference for two methods - single polarity (original) and dual
% polarities (optimized) GIRF calculation.
% Author: Zhe "Tim" Wu
% Created: Nov. 1, 2021

%% User defined parameters
% Select which gradient axis for GIRF calculation
% Select from 'X', 'Y', and 'Z'.
if exist('gradientAxis', 'var') ~= 1
    gradientAxis = 'x';
end

% Set the data path that stores subfolders 'CalculatedGIRF'
if exist('dataPath', 'var') ~= 1
    dataPath = '../DataISMRM2022';
end

% Path to load the pre-calculated GIRF
preCalcGIRFPath = strcat(dataPath, '/CalculatedGIRF/');

%% File name and path
gradientAxis = lower(gradientAxis);

resultFileName1 = strcat('GIRFOrigin_G', gradientAxis, '_Meas2.mat');
resultFileName2 = strcat('GIRFOptimized_G', gradientAxis, '_Meas2.mat');

% This will load the following variables:
% GIRF_FT, dwellTimeSig, isAvgRepetition, roPts, roTime
load(strcat(preCalcGIRFPath, resultFileName1));
GIRF_FT1 = GIRF_FT; clear GIRF_FT;
load(strcat(preCalcGIRFPath, resultFileName2));
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

maxSNRRatio = max(SNRSmooth2) ./ max(SNRSmooth1);
disp(['Ratio between max SNRs on G', gradientAxis, ' is ', num2str(maxSNRRatio)]);

%% Plot SNR
freq_fullrange = 1 / (dwellTimeSig / 1e6) / 1e3; % Full spectrum width, in unit of kHz
freq = linspace(-freq_fullrange/2, freq_fullrange/2, roPts);
freq = freq(:);

dispFreqRange = [-30, 30]; % in unit of kHz

figure(555);
set(gcf,'color','white');
plot(freq, SNRSmooth1, 'r', 'LineWidth', 1);
hold on;
plot(freq, SNRSmooth2, 'k', 'LineWidth', 1);
hold on;
xlim(dispFreqRange);
xlabel('Frequency [kHz]','FontSize', 14); ylabel('SNR [AU]','FontSize', 14);
title('SNR of GIRF in Frequency Domain','FontSize', 18);
hold off;
legend('Origin', 'Optimized', 'FontSize', 14);

% Mean of SNR ratio between display range
% Noisy part < -30 kHz or > 30 kHz was removed
[~, freqDispIndexStart] = min(abs(freq - dispFreqRange(1)));
[~, freqDispIndexEnd] = min(abs(freq - dispFreqRange(end)));
freqDispIndex = freqDispIndexStart : freqDispIndexEnd;
freqDispIndex = freqDispIndex(:);
meanSNRRatio = mean(SNRSmooth2(freqDispIndex) ./ SNRSmooth1(freqDispIndex));
disp(['Ratio between mean SNRs on G', gradientAxis, ' is ', num2str(meanSNRRatio)]);
