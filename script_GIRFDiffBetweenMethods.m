% script_GIRFDiffBetweenMethods.m
% Script comparing calculated GIRF with the original (positive blip polarity) and optimized (dual blip polarity) methods.
% This script generates subfigures on the first row in Figure 4

% Author: Zhe "Tim" Wu
% Created: Nov 1, 2021
%
%    Copyright (C) 2021-2022
%    Brain Research in Advanced Imaging and Neuromodeling - Toronto (BRAIN-TO) Lab
%    Techna Institute
%    University Health Network
%    Please see LICENSE file for details on usage

%% User defined parameters
% Select which gradient axis for GIRF calculation
% Select from 'X', 'Y', and 'Z'.
if exist('gradientAxis', 'var') ~= 1
    gradientAxis = 'z';
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
figure(444);
clf;
set(gcf,'color','white');
plot(freq, abs(GIRF_FT_mean1), 'r', 'LineWidth', 2);
hold on;
plot(freq, abs(GIRF_FT_mean2), 'k', 'LineWidth', 2);
hold on;
plot(freq, abs(GIRF_FT_mean1 - GIRF_FT_mean2), 'b', 'LineWidth', 1);
hold on;
xlim(dispFreqRange); ylim([0, 1.1]);
xlabel('Frequency [kHz]','FontSize', 18); ylabel('Magnitude of GIRF','FontSize', 18);
title(strcat('GIRF in Frequency Domain for G', gradientAxis), 'FontSize', 22);
legend('Origin', 'Optimized', 'Difference','FontSize', 18);
