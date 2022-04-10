% Main script for calculating GIRF with phantom-based method from T2*

%% User Settings: data path and gradient axis

% Set the data path that stores subfolders 'meas1' and 'meas2'
dataPath = '../DataISMRM2022';

% Select which gradient axis for GIRF calculation
% Select from 'X', 'Y', and 'Z' (case).
gradientAxis = 'z';

%% Check the validities of user inputs
if ~strcmp(gradientAxis, 'X') && ~strcmp(gradientAxis, 'x') && ...
        ~strcmp(gradientAxis, 'Y') && ~strcmp(gradientAxis, 'y') && ...
        ~strcmp(gradientAxis, 'Z') && ~strcmp(gradientAxis, 'z')
    error('Selected gradient axis must be X/x or Y/y or Z/z.');
end

gradientAxis = lower(gradientAxis);

%% Include path for utility functions
addpath('./utils/');

%% Setting data path and file names according to user selections
% Only the 2nd measurement contains the optimized protocol of GIRF.
fullDataPath = strcat(dataPath, '/Meas2/');

% We need four files for each gradient axis,
% incuding raw T2* decay signal files and ref scans for both slices.
% For Gx axis, use files with L and R as slice 1 and 2, respectively.
% For Gy axis, use files with P and A as slice 1 and 2, respectively.
% For Gz axis, use files with H and F as slice 1 and 2, respectively.
% Slices are all in 1mm thickness and 34 mm slice gap. ADC dwell time 5us.

switch gradientAxis
    % Slice Pair 1: For Gx GIRF, we need two slices along the sagittal plane.
    case {'X', 'x'}
        fn_Slice1_POS = 'PositiveSagL17.mat';
        fn_Slice1_NEG = 'NegativeSagL17.mat';
        fn_Slice2_POS = 'PositiveSagR17.mat';
        fn_Slice2_NEG = 'NegativeSagR17.mat';
        fn_Slice1_ref = 'RefSagL17.mat';
        fn_Slice2_ref = 'RefSagR17.mat';
    
    % Slice Pair 2: For Gy GIRF, we need two slices along the coronal plane.
    case {'Y', 'y'}
        fn_Slice1_POS = 'PositiveCorP17.mat';
        fn_Slice1_NEG = 'NegativeCorP17.mat';
        fn_Slice2_POS = 'PositiveCorA17.mat';
        fn_Slice2_NEG = 'NegativeCorA17.mat';
        fn_Slice1_ref = 'RefCorP17.mat';
        fn_Slice2_ref = 'RefCorA17.mat';

    % Slice Pair 3: For Gz GIRF, we need two slices along the transversal plane.
    case {'Z', 'z'}
        fn_Slice1_POS = 'PositiveTraH17.mat';
        fn_Slice1_NEG = 'NegativeTraH17.mat';
        fn_Slice2_POS = 'PositiveTraF17.mat';
        fn_Slice2_NEG = 'NegativeTraF17.mat';
        fn_Slice1_ref = 'RefTraH17.mat';
        fn_Slice2_ref = 'RefTraF17.mat';
        
    otherwise
        error('Selected gradient axis must be X/x or Y/y or Z/z.');
end

% File containing all gradient waveforms (blips)
fn_gradient = 'InputGradients.mat';

fn_Slice1_POS = strcat(fullDataPath, fn_Slice1_POS);
fn_Slice1_NEG = strcat(fullDataPath, fn_Slice1_NEG);
fn_Slice2_POS = strcat(fullDataPath, fn_Slice2_POS);
fn_Slice2_NEG = strcat(fullDataPath, fn_Slice2_NEG);
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
[rawSigS1_POS, rawSigS1_NEG, rawSigS2_POS, rawSigS2_NEG, refS1, refS2] = readMultiFiles('kspace_all', fn_Slice1_POS, fn_Slice1_NEG, fn_Slice2_POS, fn_Slice2_NEG, fn_Slice1_ref, fn_Slice2_ref);
% Load all gradient blips, includes gradIn_all
load(fn_gradient, 'gradIn_all'); 

%% Step 2: Processing gradient inputs
params.roPts = size(rawSigS1_POS, 1); % Readout points
params.nRep = size(rawSigS1_POS, 2); % Number of repetition
params.nGradAmp = size(rawSigS1_POS, 3); % Number of gradient blips

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
rawSigS1_POS = permute(rawSigS1_POS, [1, 3, 2]); % [nRO, nGradAmp, nRep]
rawSigS1_NEG = permute(rawSigS1_NEG, [1, 3, 2]);
rawSigS2_POS = permute(rawSigS2_POS, [1, 3, 2]);
rawSigS2_NEG = permute(rawSigS2_NEG, [1, 3, 2]);

% Average the ref scan data across the average dimension
refS1 = squeeze(mean(refS1,2)); % [nRO, nRefScan]
refS2 = squeeze(mean(refS2,2));

[gradOutput, gradOutputFT] = calculateOutputGradientOptimized(rawSigS1_POS, rawSigS1_NEG, rawSigS2_POS, rawSigS2_NEG, refS1, refS2, params);

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

displayGIRFMagnitude(GIRF_FT, freqFull, dispFreqRange, 555);

%% Step 5: Save Results
if exist(dataSavePath, 'file') ~= 7
    mkdir(dataSavePath);
end

save(strcat(dataSavePath, 'GIRFOptimized_G', gradientAxis, '_Meas2'), 'GIRF_FT', 'params', 'roTime');
