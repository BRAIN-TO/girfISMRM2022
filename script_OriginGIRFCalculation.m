% Main script for calculating GIRF with phantom-based method from T2*

%% Include path for utility functions
addpath('./utils/');

%% Read Measured Data

dataPath = '../DataISMRM2022/Meas1/';
% dataPath = '../DataISMRM2022/Meas2/';

% We need four files for each gradient axis,
% incuding raw T2* decay signal files and ref scans for both slices.
% For Gx axis, use files with L and R as slice 1 and 2, respectively.
% For Gy axis, use files with P and A as slice 1 and 2, respectively.
% For Gz axis, use files with H and F as slice 1 and 2, respectively.
% Slices are all in 1mm thickness and 34 mm slice gap. ADC dwell time 5us.

% % Slice Pair 1: Sagittal slices (Along X) 
% fn_Slice1 = 'PosSignalSagL17.mat';
% fn_Slice2 = 'PosSignalSagR17.mat';
% fn_Slice1_ref = 'RefSagL17.mat';
% fn_Slice2_ref = 'RefSagR17.mat';

% % Slice Pair 2: Coronal slices (Along Y)
% fn_Slice1 = 'PosSignalCorP17.mat';
% fn_Slice2 = 'PosSignalCorA17.mat';
% fn_Slice1_ref = 'RefCorP17.mat';
% fn_Slice2_ref = 'RefCorA17.mat';

% Slice Pair 3: Transverse slices (Along Z)
fn_Slice1 = 'PosSignalTraH17.mat';
fn_Slice2 = 'PosSignalTraF17.mat';
fn_Slice1_ref = 'RefTraH17.mat';
fn_Slice2_ref = 'RefTraF17.mat';

% File containing all gradient waveforms (blips)
fn_gradient = 'gradients_zerofilled.mat';

fn_Slice1 = strcat(dataPath, fn_Slice1);
fn_Slice2 = strcat(dataPath, fn_Slice2);
fn_Slice1_ref = strcat(dataPath, fn_Slice1_ref);
fn_Slice2_ref = strcat(dataPath, fn_Slice2_ref);
fn_gradient = strcat(dataPath, fn_gradient);

%% Fixed parameters
params = [];
params.gammabar = 42.57e3; % in unit of Hz/mT
params.slicePos = 0.017; % distance of both slices from isocenter, in unit of m
params.gradRasterTime = 10; % in us
params.adcDwellTime = 5; % in us

%% Start Data processing

%% Load all necessary files

% Load T2* decay signal
[rawSigS1, rawSigS2, refS1, refS2] = readMultiFiles(fn_Slice1, fn_Slice2, fn_Slice1_ref, fn_Slice2_ref);
% Load all gradient blips, includes gradIn_all
load(fn_gradient, 'gradIn_all'); 

%% Processing gradient inputs
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


%% Calculate Output Gradients and GIRF

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

%% Plot GIRF in Frequency Domain

% Calculate full frequency range of the GIRF spectrum
freqRange = 1 / (params.adcDwellTime / 1e6) / 1e3; % Full spectrum width, in unit of kHz
freqFull = linspace(-freqRange/2, freqRange/2, params.roPts);
freqFull = freqFull(:);

dispFreqRange = [-30, 30]; % in unit of kHz

displayGIRFMagnitude(GIRF_FT, freqFull, dispFreqRange, 555);
