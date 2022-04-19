% script_OriginGIRFCalculation.m
% Script for calculate GIRF with phantom-based method from T2* using original method (positive triangle blips only)
% This script is using the coil-averaged data for demonstration purpose ONLY.
% Please use the pre-calculated GIRF from the non-compressed 

% Author: Zhe "Tim" Wu
% Created: Nov 1, 2021
%
%    Copyright (C) 2021-2022
%    Brain Research in Advanced Imaging and Neuromodeling - Toronto (BRAIN-TO) Lab
%    Techna Institute
%    University Health Network
%    Please see LICENSE file for details on usage

%% User Settings: data path and gradient axis

% Set the data path that stores subfolders 'meas1' and 'meas2'
if exist('dataPath', 'var') ~= 1
    dataPath = '../DataISMRM2022';
end

% Set the path for saving GIRF results
dataSavePath = strcat(dataPath, '/Results/');

% We have two measurements for stability investigations with a 8-month gap
% Select from 1 and 2.
if exist('measNum', 'var') ~= 1
    measNum = 1;
end

% Select which gradient axis for GIRF calculation
% Select from 'X', 'Y', and 'Z'.
if exist('gradientAxis', 'var') ~= 1
    gradientAxis = 'x';
end

%% Check the validities of user inputs
if measNum ~= 1 && measNum ~= 2
    error('Measurement number (measNum) must be 1 or 2.');
end

if ~strcmp(gradientAxis, 'X') && ~strcmp(gradientAxis, 'x') && ...
        ~strcmp(gradientAxis, 'Y') && ~strcmp(gradientAxis, 'y') && ...
        ~strcmp(gradientAxis, 'Z') && ~strcmp(gradientAxis, 'z')
    error('Selected gradient axis must be X/x or Y/y or Z/z.');
end

gradientAxis = lower(gradientAxis);

%% Include path for utility functions
addpath('./utils/');

%% Setting data path and file names according to user selections
fullDataPath = strcat(dataPath, '/Meas', num2str(measNum), '/');

% We need four files for each gradient axis,
% incuding raw T2* decay signal files and ref scans for both slices.
% For Gx axis, use files with L and R as slice 1 and 2, respectively.
% For Gy axis, use files with P and A as slice 1 and 2, respectively.
% For Gz axis, use files with H and F as slice 1 and 2, respectively.
% Slices are all in 1mm thickness and 34 mm slice gap. ADC dwell time 5us.

switch gradientAxis
    % Slice Pair 1: For Gx GIRF, we need two slices along the sagittal plane.
    case {'X', 'x'}
        fn_Slice1 = 'PositiveSagL17.mat';
        fn_Slice2 = 'PositiveSagR17.mat';
        fn_Slice1_ref = 'RefSagL17.mat';
        fn_Slice2_ref = 'RefSagR17.mat';
    
    % Slice Pair 2: For Gy GIRF, we need two slices along the coronal plane.
    case {'Y', 'y'}
        fn_Slice1 = 'PositiveCorP17.mat';
        fn_Slice2 = 'PositiveCorA17.mat';
        fn_Slice1_ref = 'RefCorP17.mat';
        fn_Slice2_ref = 'RefCorA17.mat';

    % Slice Pair 3: For Gz GIRF, we need two slices along the transversal plane.
    case {'Z', 'z'}
        fn_Slice1 = 'PositiveTraH17.mat';
        fn_Slice2 = 'PositiveTraF17.mat';
        fn_Slice1_ref = 'RefTraH17.mat';
        fn_Slice2_ref = 'RefTraF17.mat';
    
    otherwise
        error('Selected gradient axis must be X/x or Y/y or Z/z.');
end

% File containing all gradient waveforms (blips)
fn_gradient = 'InputGradients.mat';

fn_Slice1 = strcat(fullDataPath, fn_Slice1);
fn_Slice2 = strcat(fullDataPath, fn_Slice2);
fn_Slice1_ref = strcat(fullDataPath, fn_Slice1_ref);
fn_Slice2_ref = strcat(fullDataPath, fn_Slice2_ref);
fn_gradient = strcat(fullDataPath, fn_gradient);

%% Fixed parameters
params = [];
params.gammabar = 42.57e3; % in unit of Hz/mT
params.slicePos = 0.017; % distance of both slices from isocenter, in unit of m
params.gradRasterTime = 10; % in us
params.adcDwellTime = 5; % in us

%% Start Data processing

%% Step 1: Load all necessary files

% Load T2* decay signal. The variable name in each MAT file is called 'kspace_all'
[rawSigS1, rawSigS2, refS1, refS2] = readMultiFiles('kspace_all', fn_Slice1, fn_Slice2, fn_Slice1_ref, fn_Slice2_ref);
% Load all gradient blips, includes gradIn_all
load(fn_gradient, 'gradIn_all'); 

%% Step 2: Processing gradient inputs
params.roPts = size(rawSigS1, 1); % Readout points
params.nRep = size(rawSigS1, 2); % Number of repetition
params.nGradAmp = size(rawSigS1, 3); % Number of gradient blips

[gResamp, roTime] = resampGradients(gradIn_all, params);

% The nominal starting time of the blips is 2000us after ADC starts
timeShift = 2000; % in us
gradInput = circshift(gResamp, timeShift/params.adcDwellTime, 1);

% Note that we have NOT expanded the input grad to nCoil dim
% We will decide this issue later during output signal processing.
gradInputFT = fftshift(fft(fftshift(gradInput,1),[],1),1); % [nRO, nGradAmp]
gradInputFT = repmat(gradInputFT, [1,1,params.nRep]);


%% Step 3: Calculate Output Gradients and GIRF

% Change dimention from [nRO, nRep, nGradAmp] to [nRO, nGradAmp, nRep]
rawSigS1 = permute(rawSigS1, [1, 3, 2]); % [nRO, nGradAmp, nRep]
rawSigS2 = permute(rawSigS2, [1, 3, 2]);

% Average the ref scan data across the average dimension
refS1 = squeeze(mean(refS1,2)); % [nRO, nRefScan]
refS2 = squeeze(mean(refS2,2));

[gradOutput, gradOutputFT] = calculateOutputGradient(rawSigS1, rawSigS2, refS1, refS2, params);

% Sum up along the last dimension, which is always the nGradAmp dimension,
% no matter whether the nCoil dimension exists or not.
GIRF_FT = sum(gradOutputFT .* conj(gradInputFT), 2) ./ sum(abs(gradInputFT).^2, 2);
GIRF_FT = squeeze(GIRF_FT);

%% Step 4: Plot GIRF in Frequency Domain

% Calculate full frequency range of the GIRF spectrum
freqRange = 1 / (params.adcDwellTime / 1e6) / 1e3; % Full spectrum width, in unit of kHz
freqFull = linspace(-freqRange/2, freqRange/2, params.roPts);
freqFull = freqFull(:);

dispFreqRange = [-30, 30]; % in unit of kHz

displayGIRFMagnitude(GIRF_FT, freqFull, dispFreqRange, 111);

%% Step 5: Save Results
if exist(dataSavePath, 'file') ~= 7
    mkdir(dataSavePath);
end

save(strcat(dataSavePath, 'GIRFOrigin_G', gradientAxis, '_Meas', num2str(measNum)), 'GIRF_FT', 'params', 'roTime');
