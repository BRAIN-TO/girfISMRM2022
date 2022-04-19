% Script main.m
% Run all the demonstration scripts for Gradient Impulse Response Function (GIRF) calculation and analysis (stability, SNR...)

% Author: Zhe "Tim" Wu
% Created: Apr 10, 2022
%
%    Copyright (C) 2021-2022
%    Brain Research in Advanced Imaging and Neuromodeling - Toronto (BRAIN-TO) Lab
%    Techna Institute
%    University Health Network
%    Please see LICENSE file for details on usage

%% Set User Parameters
clear;
clc;

% Set the data path that stores the downloaded data
dataPath = '../DataISMRM2022/';

% Select which gradient axis for GIRF calculation
% Select from 'X', 'Y', and 'Z'.
gradientAxis = 'Z';

% Only applied to script_OriginGIRFCalculation.m
% We have two measurements for stability investigations with a 8-month gap
% Select from 1 and 2.
measNum = 1;

%% Print General Information

disp(['Running All Scripts for Gradient Axis ', upper(gradientAxis), ':']);
disp(' ');

%% Demo 1 Original GIRF Calculation (Using Single Blips)
disp('Demo 1: Calculating GIRF using the original method (positive blips). Results are in Figure 111.');
disp(strcat('The data from measurement ', num2str(measNum), ' is used.'));
script_OriginGIRFCalculation;
disp(' ');

%% Demo 2 Optimized GIRF Calculation (Using Dual Blips)
disp('Demo 2: Calculating GIRF using the optimized method (both positive and negative blips). Results are in Figure 222.');
script_OptimizedGIRFCalculation;
disp(' ');

%% Demo 3 GIRF stability analysis for 8-month gap (Figure 3 C-E)
disp('Demo 3: Analysing stability of GIRF for both measurements with 8-month gap using the original method. Results are in Figure 333.');
script_GIRFStability;
disp(' ');

%% Demo 4 Comparison between two GIRF calculation methods (Figure 4 first row)
disp('Demo 4: Comparing calculated GIRF from both original and optimized method from Measurement 2. Results are in Figure 444.');
script_GIRFDiffBetweenMethods;
disp(' ');

%% Demo 5 Analyzing the SNR difference for two methods - single polarity (original) and dual polarities (optimized) GIRF calculation.
disp('Demo 5: Analyzing the SNR difference of the calculated GIRF from both original and optimized method from Measurement 2. Results are in Figure 555.');
script_SNRAnalysis;
disp(' ');
