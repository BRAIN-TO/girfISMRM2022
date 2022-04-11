% script_GIRFStability.m
%Script for GIRF stability analysis
% This script generates Figure 3 (C-E)
% The relative difference of two measurements is based on the frequency
% ranges with 99% energy concentration for spiral and EPI trajectories. See
% Figure 3 (A-B)
% Author: Zhe "Tim" Wu
% Created: Nov 1, 2021

%% User defined parameters
% Select which gradient axis for GIRF calculation
% Select from 'X', 'Y', and 'Z'.
gradientAxis = 'x';

% Path to load the pre-calculated GIRF
resultPath = "../DataISMRM2022/CalculatedGIRF/";

%% File name and path
gradientAxis = lower(gradientAxis);

resultFileName1 = strcat('GIRFOrigin_G', gradientAxis, '_Meas1.mat');
resultFileName2 = strcat('GIRFOrigin_G', gradientAxis, '_Meas2.mat');

% This will load the following variables:
% GIRF_FT, dwellTimeSig, isAvgRepetition, roPts, roTime
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

GIRF_FT_mean1 = abs(mean(GIRF_FT1,2));
GIRF_FT_mean2 = abs(mean(GIRF_FT2,2));

% Frequency range to observe in kHz (contains 99% of the energy of spiral and EPI trajectories)
% See Figure 3 (A-B)
freq1 = 3.2; 
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
xlim(dispFreqRange); ylim([0, 1.1]);
xlabel('Frequency [kHz]','FontSize', 18); ylabel('Magnitude of GIRF','FontSize', 18);
title(strcat('GIRF of G', gradientAxis, ' in the Frequency Domain'), 'FontSize', 22);
legend('Initial', '8 Months Later', 'Difference','FontSize', 18);

girfRelDiff = abs(GIRF_FT_mean1 - GIRF_FT_mean2) ./ GIRF_FT_mean1;
girfRelDiffFreq1 = max(girfRelDiff(indexLeft1:indexRight1));
girfRelDiffFreq2 = max(girfRelDiff(indexLeft2:indexRight2));

disp(['The max relative difference of GIRF on G', gradientAxis, ' in the frequency range of ±', num2str(freq1), 'kHz is ', num2str(girfRelDiffFreq1 * 100), '%']);
disp(['The max relative difference of GIRF on G', gradientAxis, ' in the frequency range of ±', num2str(freq2), 'kHz is ', num2str(girfRelDiffFreq2 * 100), '%']);
