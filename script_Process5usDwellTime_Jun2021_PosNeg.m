%% Save File Path
resultPath = "../ISMRM2022Results/";
resultFileName = "2021Jun_PosNeg_Gx.mat";

%% Read Measured Data
% % The dataset for Gx Gy and Gz, 1mm thickness, 34 mm slice gap, 5 us dwell time
fpath_mat = '../Data20210629/'; 

% % Slice Pair 1: Transverse slices (Along Z)
% fn_mat_Slice1_POS = 'PosSignalTraH17_5us.mat';
% fn_mat_Slice2_POS = 'PosSignalTraF17_5us.mat';
% fn_mat_Slice1_NEG = 'NegSignalTraH17_5us.mat';
% fn_mat_Slice2_NEG = 'NegSignalTraF17_5us.mat';
% fn_mat_Slice1_ref = 'RefTraH17_5us.mat';
% fn_mat_Slice2_ref = 'RefTraF17_5us.mat';

% % Slice Pair 2: Coronal slices (Along Y)
% fn_mat_Slice1_POS = 'PosSignalCorP17_5us.mat';
% fn_mat_Slice2_POS = 'PosSignalCorA17_5us.mat';
% fn_mat_Slice1_NEG = 'NegSignalCorP17_5us.mat';
% fn_mat_Slice2_NEG = 'NegSignalCorA17_5us.mat';
% fn_mat_Slice1_ref = 'RefCorP17_5us.mat';
% fn_mat_Slice2_ref = 'RefCorA17_5us.mat';

% Slice Pair 3: Sagittal slices (Along X)
fn_mat_Slice1_POS = 'PosSignalSagL17_5us.mat';
fn_mat_Slice2_POS = 'PosSignalSagR17_5us.mat';
fn_mat_Slice1_NEG = 'NegSignalSagL17_5us.mat';
fn_mat_Slice2_NEG = 'NegSignalSagR17_5us.mat';
fn_mat_Slice1_ref = 'RefSagL17_5us.mat';
fn_mat_Slice2_ref = 'RefSagR17_5us.mat';

% Containing variables: 'kspace_all', 'roPts', 'nch', 'avgNum', 'acqNum', 'gradAmp', 'dwellTime', 'roTime'
load(strcat(fpath_mat, fn_mat_Slice1_POS));
rawSigS1_POS = kspace_all; clear kspace_all;
load(strcat(fpath_mat, fn_mat_Slice1_NEG));
rawSigS1_NEG = kspace_all; clear kspace_all;

% Containing variables: 'kspace_all', 'roPts', 'nch', 'avgNum', 'acqNum', 'gradAmp', 'dwellTime', 'roTime'
load(strcat(fpath_mat, fn_mat_Slice2_POS));
rawSigS2_POS = kspace_all; clear kspace_all;
load(strcat(fpath_mat, fn_mat_Slice2_NEG));
rawSigS2_NEG = kspace_all; clear kspace_all;

% Containing variables: 'kspace_all', 'roPts', 'nch', 'avgNum', 'acqNum', 'gradAmp', 'dwellTime', 'roTime'
load(strcat(fpath_mat, fn_mat_Slice1_ref));
refS1 = kspace_all; clear kspace_all;

% Containing variables: 'kspace_all', 'roPts', 'nch', 'avgNum', 'acqNum', 'gradAmp', 'dwellTime', 'roTime'
load(strcat(fpath_mat, fn_mat_Slice2_ref));
refS2 = kspace_all; clear kspace_all;

% Confirm the size of data
disp(['Size of raw signal of Slice 1 Positive Polarity: ', num2str(size(rawSigS1_POS))]);
disp(['Size of raw signal of Slice 1 Negative Polarity: ', num2str(size(rawSigS1_NEG))]);
disp(['Size of raw signal of Slice 2 Positive Polarity: ', num2str(size(rawSigS2_POS))]);
disp(['Size of raw signal of Slice 2 Negative Polarity: ', num2str(size(rawSigS2_NEG))]);
disp(['Size of ref data of Slice 1: ', num2str(size(refS1))]);
disp(['Size of ref data of Slice 2: ', num2str(size(refS2))]);

%% Truncate the ADC window to match Nov 2020 Data
rawSigS1_POS = rawSigS1_POS(1:10002, :, :, :);
rawSigS1_NEG = rawSigS1_NEG(1:10002, :, :, :);
rawSigS2_POS = rawSigS2_POS(1:10002, :, :, :);
rawSigS2_NEG = rawSigS2_NEG(1:10002, :, :, :);
refS1 = refS1(1:10002, :, :, :);
refS2 = refS2(1:10002, :, :, :);
roPts = 10002;
roTime = roTime(1:10002);

%% Read Input Gradient Files
fpath_mat = '../GIRF_Gradient/mat/';
fn_mat = 'gradients_zerofilled.mat';

load(strcat(fpath_mat, fn_mat)); % gradIn_all, nSamplesPerAxis

%% Data processing

%% Processing gradient input
dwellTimeGrad = 10; % in us
dwellTimeSig = 5; % in us

% Interpolate the input gradient and upsample it from dwellTime = 10us to 5us
gradTime = 0:(nSamplesPerAxis - 1);
gradTime = gradTime * dwellTimeGrad;
gradTime = gradTime(:);
roTime = roTime(:);

nGradIn = size(gradIn_all,2); % Number of input gradients

% We have already have the timeline for raw signal as roTime
gradIn_all_interp = [];
for fc = 1 : nGradIn
    gradIn_all_interp = cat(2, gradIn_all_interp, interp1(gradTime, gradIn_all(:,fc), roTime));
end
gradIn_all_interp(end,:) = 0; % Fix the interp error at the last point.

% Add zeros at the beginning of gradient for 2000us, and truncate the last
% 2000 us
timeShift = 2000; % in us
gradInput = circshift(gradIn_all_interp, timeShift/dwellTimeSig, 1);

% Note that we have NOT expanded the input grad to nCoil dim
% We will decide this issue later during output signal processing.
gradInputFT = fftshift(fft(fftshift(gradInput,1),[],1),1); % [nRO, nGradAmp]

%% Calculate GIRF with Method 1: first temporal derivitive of phase, then subtracte phases from two slices

isAvgRepetition = 0; % Whether we use mean across dim of repetitions of raw signal
isAvgAcrossCoil = 1; % Whether we use mean across dim nCoil or preserve all the coils.
useRefData = 1; % Whether we are using reference T2* decay to correct phase drifting

gammabar = 42.57e3; % in unit of Hz/mT
slicePos = 0.017; % in unit of m

% Average the raw signal across the repetition dimension
if isAvgRepetition
    avgSigS1_POS = squeeze(mean(rawSigS1_POS,3)); % [nRO, nCoil, nGradAmp]
    avgSigS1_NEG = squeeze(mean(rawSigS1_NEG,3)); % [nRO, nCoil, nGradAmp]
    avgSigS2_POS = squeeze(mean(rawSigS2_POS,3)); % [nRO, nCoil, nGradAmp]
    avgSigS2_NEG = squeeze(mean(rawSigS2_NEG,3)); % [nRO, nCoil, nGradAmp]
else
    avgSigS1_POS = permute(rawSigS1_POS, [1, 2, 4, 3]); % [nRO, nCoil, nGradAmp, nRep]
    avgSigS1_NEG = permute(rawSigS1_NEG, [1, 2, 4, 3]); % [nRO, nCoil, nGradAmp, nRep]
    avgSigS2_POS = permute(rawSigS2_POS, [1, 2, 4, 3]); % [nRO, nCoil, nGradAmp, nRep]
    avgSigS2_NEG = permute(rawSigS2_NEG, [1, 2, 4, 3]); % [nRO, nCoil, nGradAmp, nRep]
end

if useRefData
    % Average the ref scan data across the average dimension
    % Remember: We have ref scan before every two blip measurements
    % An additional one ref scan was added to the end of the whole
    % So nRefScan = nGradAmp / 2 + 1
    avgRefS1 = squeeze(mean(refS1,3)); % [nRO, nCoil, nRefScan]
    avgRefS2 = squeeze(mean(refS2,3));
    
    % Need to repeate if raw signal is not averaged
    % Note we are using mean of ref signal because they are independent to
    % the measurement of raw signal decay with blips; besides, the
    % repetition numbers are different in ref signal acquisition (5, vs 50 in raw signal)
    nRep = size(avgSigS1_POS,4);
    if isAvgRepetition == 0
        avgRefS1 = repmat(avgRefS1, [1,1,1,nRep]);
        avgRefS2 = repmat(avgRefS2, [1,1,1,nRep]);
    end

    avgSigS1Corr_POS = zeros(size(avgSigS1_POS));
    avgSigS1Corr_NEG = zeros(size(avgSigS1_NEG));
    avgSigS2Corr_POS = zeros(size(avgSigS2_POS));
    avgSigS2Corr_NEG = zeros(size(avgSigS2_NEG));

    nGradAmp = length(gradAmp);

    for nn = 1 : nGradAmp
        if isAvgRepetition
            avgSigS1Corr_POS(:,:,nn) = avgSigS1_POS(:,:,nn) ./ avgRefS1(:,:,floor((nn-1)/3)+1);
            avgSigS1Corr_NEG(:,:,nn) = avgSigS1_NEG(:,:,nn) ./ avgRefS1(:,:,floor((nn-1)/3)+1);
            avgSigS2Corr_POS(:,:,nn) = avgSigS2_POS(:,:,nn) ./ avgRefS2(:,:,floor((nn-1)/3)+1);
            avgSigS2Corr_NEG(:,:,nn) = avgSigS2_NEG(:,:,nn) ./ avgRefS2(:,:,floor((nn-1)/3)+1);
        else
            avgSigS1Corr_POS(:,:,nn,:) = avgSigS1_POS(:,:,nn,:) ./ avgRefS1(:,:,floor((nn-1)/3)+1,:);
            avgSigS1Corr_NEG(:,:,nn,:) = avgSigS1_NEG(:,:,nn,:) ./ avgRefS1(:,:,floor((nn-1)/3)+1,:);
            avgSigS2Corr_POS(:,:,nn,:) = avgSigS2_POS(:,:,nn,:) ./ avgRefS2(:,:,floor((nn-1)/3)+1,:);
            avgSigS2Corr_NEG(:,:,nn,:) = avgSigS2_NEG(:,:,nn,:) ./ avgRefS2(:,:,floor((nn-1)/3)+1,:);
        end
    end

    avgSigS1Diff_POS = diff(unwrap(angle(avgSigS1Corr_POS),[],1),1,1) / (dwellTimeSig / 1e6); % in rad/s
    avgSigS1Diff_NEG = diff(unwrap(angle(avgSigS1Corr_NEG),[],1),1,1) / (dwellTimeSig / 1e6); % in rad/s
    avgSigS2Diff_POS = diff(unwrap(angle(avgSigS2Corr_POS),[],1),1,1) / (dwellTimeSig / 1e6);
    avgSigS2Diff_NEG = diff(unwrap(angle(avgSigS2Corr_NEG),[],1),1,1) / (dwellTimeSig / 1e6);
else
    avgSigS1Diff_POS = diff(unwrap(angle(avgSigS1_POS),[],1),1,1) / (dwellTimeSig / 1e6); % in rad/s
    avgSigS1Diff_NEG = diff(unwrap(angle(avgSigS1_NEG),[],1),1,1) / (dwellTimeSig / 1e6); % in rad/s
    avgSigS2Diff_POS = diff(unwrap(angle(avgSigS2_POS),[],1),1,1) / (dwellTimeSig / 1e6);
    avgSigS2Diff_NEG = diff(unwrap(angle(avgSigS2_NEG),[],1),1,1) / (dwellTimeSig / 1e6);
end

if isAvgRepetition
    avgSigS1Diff_POS = cat(1, zeros(1, size(avgSigS1Diff_POS,2), size(avgSigS1Diff_POS,3)), avgSigS1Diff_POS);
    avgSigS1Diff_NEG = cat(1, zeros(1, size(avgSigS1Diff_NEG,2), size(avgSigS1Diff_NEG,3)), avgSigS1Diff_NEG);
    avgSigS2Diff_POS = cat(1, zeros(1, size(avgSigS2Diff_POS,2), size(avgSigS2Diff_POS,3)), avgSigS2Diff_POS);
    avgSigS2Diff_NEG = cat(1, zeros(1, size(avgSigS2Diff_NEG,2), size(avgSigS2Diff_NEG,3)), avgSigS2Diff_NEG);
else
    avgSigS1Diff_POS = cat(1, zeros(1, size(avgSigS1Diff_POS,2), size(avgSigS1Diff_POS,3), size(avgSigS1Diff_POS,4)), avgSigS1Diff_POS);
    avgSigS1Diff_NEG = cat(1, zeros(1, size(avgSigS1Diff_NEG,2), size(avgSigS1Diff_NEG,3), size(avgSigS1Diff_NEG,4)), avgSigS1Diff_NEG);
    avgSigS2Diff_POS = cat(1, zeros(1, size(avgSigS2Diff_POS,2), size(avgSigS2Diff_POS,3), size(avgSigS2Diff_POS,4)), avgSigS2Diff_POS);
    avgSigS2Diff_NEG = cat(1, zeros(1, size(avgSigS2Diff_NEG,2), size(avgSigS2Diff_NEG,3), size(avgSigS2Diff_NEG,4)), avgSigS2Diff_NEG);
end

gradOutput = ((avgSigS1Diff_POS - avgSigS1Diff_NEG) + (avgSigS2Diff_NEG - avgSigS2Diff_POS)) / 4; % Phase averaged from two slices, in unit of rad/s
gradOutput = gradOutput / (gammabar * 2*pi) / slicePos; % in unit of mT/m

b0ecOutput = ((avgSigS1Diff_POS - avgSigS1Diff_NEG) - (avgSigS2Diff_NEG - avgSigS2Diff_POS)) / 4;

if isAvgAcrossCoil == 1
    gradInputFT = fftshift(fft(fftshift(gradInput,1),[],1),1); % [nRO, nGradAmp]
    if isAvgRepetition == 0
        gradInputFT = repmat(gradInputFT, [1,1,nRep]);
    end
    gradOutput = squeeze(mean(gradOutput,2)); % Average across dim of nCoil, now matrix size is [nRO, nGradAmp]
    b0ecOutput = squeeze(mean(b0ecOutput,2)); % Average across dim of nCoil, now matrix size is [nRO, nGradAmp]
else
    gradInputFT = fftshift(fft(fftshift(gradInput,1),[],1),1); % [nRO, nGradAmp]
    % If we preserve all coils, we need to expand the input gradient along nCoil dimension
    gradInputFT = repmat(gradInputFT, [1,1,nch]); % [nRO, nGradAmp, nCoil]
    gradInputFT = permute(gradInputFT, [1,3,2]); % [nRO, nCoil, nGradAmp]
    if isAvgRepetition == 0
        gradInputFT = repmat(gradInputFT, [1,1,1,nRep]);
    end
end

gradOutputFT = fftshift(fft(fftshift(gradOutput,1),[],1),1);
b0ecOutputFT = fftshift(fft(fftshift(b0ecOutput,1),[],1),1);

% Sum up along the last dimension, which is always the nGradAmp dimension,
% no matter whether the nCoil dimension exists or not.
if isAvgRepetition
    GIRF_FT = sum(gradOutputFT .* conj(gradInputFT), ndims(gradOutputFT)) ./ sum(abs(gradInputFT).^2, ndims(gradOutputFT));
    b0ec_FT = sum(b0ecOutputFT .* conj(gradInputFT), ndims(b0ecOutputFT)) ./ sum(abs(gradInputFT).^2, ndims(b0ecOutputFT));
else
    GIRF_FT = sum(gradOutputFT .* conj(gradInputFT), 2) ./ sum(abs(gradInputFT).^2, 2);
    GIRF_FT = squeeze(GIRF_FT);
    b0ec_FT = sum(b0ecOutputFT .* conj(gradInputFT), 2) ./ sum(abs(gradInputFT).^2, 2);
    b0ec_FT = squeeze(b0ec_FT);
end

%% Save Results
save(strcat(resultPath, resultFileName), 'GIRF_FT', 'b0ec_FT', 'dwellTimeSig', 'roPts', 'isAvgRepetition', 'roTime');

%% Filtering GIRF_FT
%% (1) Truncate and Zero-filling
load(strcat(resultPath, resultFileName));

if isAvgRepetition
    GIRF_FT2 = GIRF_FT; % GIRF_FT has been averaged
else
    GIRF_FT2 = mean(GIRF_FT,2); % GIRF_FT has not been averaged so we need to do it here.
end

GIRF_FT2(1:2000,:,:) = 0;
GIRF_FT2(8001:end,:,:) = 0;
GIRF_Temporal = ifftshift(ifft(ifftshift(GIRF_FT2,1),[],1),1);

%% (2) Truncate, Zero-filling, and extend the sequence to 50000 points
% GIRF_FT2 = GIRF_FT;
% GIRF_FT2(1:2000,:,:) = 0;
% GIRF_FT2(8001:end,:,:) = 0;
% GIRF_FT3 = zeros(50000,1);
% GIRF_FT3(20001:30002) = GIRF_FT2;
% GIRF_Temporal = ifftshift(ifft(ifftshift(GIRF_FT3,1),[],1),1);

%% (3) Hann filter
% w = hann(size(GIRF_FT,1), 'periodic');
% GIRF_FT2 = GIRF_FT .* repmat(w, [1, size(GIRF_FT,2), size(GIRF_FT,3)]);
% GIRF_Temporal = ifftshift(ifft(ifftshift(GIRF_FT2,1),[],1),1);

%% Plot GIRF in Frequency Domain
freq_fullrange = 1 / (dwellTimeSig / 1e6) / 1e3; % Full spectrum width, in unit of kHz
freq = linspace(-freq_fullrange/2, freq_fullrange/2, roPts);
freq = freq(:);

dispFreqRange = [-30, 30]; % in unit of kHz

% freq_index = 4551:5451; % -9kHz ~ +9kHz
% freq_index = 4501:5501; % -10kHz ~ +10kHz
% freq_index = 4251:5751; % -15kHz ~ +15kHz
% freq_index = 4001:6001; % -20kHz ~ +20kHz
% freq_index = 3501:6501; % -30kHz ~ +30kHz
% freq_index = 3001:7001; % -40kHz ~ +40kHz

figure(666);
set(gcf,'color','white');

if isAvgRepetition
    subplot(211);
    plot(freq, abs(GIRF_FT));
    xlim(dispFreqRange); ylim([0, 1.1]);
    xlabel('Frequency [kHz]'); ylabel('Magnitude of GIRF [AU]');
    title('Magnitude of GIRF in Frequency Domain');

    subplot(212);
    plot(freq, angle(GIRF_FT));
    xlim(dispFreqRange); ylim([-4,4]);
    xlabel('Frequency [kHz]'); ylabel('Phase of GIRF [rad]');
    title('Phase of GIRF in Frequency Domain');
else
    [~, freqDispIndexStart] = min(abs(freq - dispFreqRange(1)));
    [~, freqDispIndexEnd] = min(abs(freq - dispFreqRange(end)));
    freqDispIndex = freqDispIndexStart : freqDispIndexEnd;
    freqDispIndex = freqDispIndex(:);
    
    GIRF_FT_mean = mean(GIRF_FT,2);
    GIRF_FT_std = std(GIRF_FT, 0, 2);
    GIRF_FT_mean_abs = abs(GIRF_FT_mean);
    GIRF_FT_mean_phase = angle(GIRF_FT_mean);
%     GIRF_FT_std_abs = std(abs(GIRF_FT), 0, 2);
%     GIRF_FT_std_phase = std(angle(GIRF_FT), 0, 2);
    GIRF_FT_std_abs = abs(GIRF_FT_std);
    GIRF_FT_std_phase = angle(GIRF_FT_std);
    
    subplot(211);
%     errorBarColor = [255 222 173]/255; % navajoWhite
    errorBarColor = [255 99 71]/255; % Tomato Red
%     errorBarColor = [128 193 219]/255; % shallow blue

    patch = fill([freq',fliplr(freq')], [GIRF_FT_mean_abs' - GIRF_FT_std_abs', fliplr(GIRF_FT_mean_abs' + GIRF_FT_std_abs')], errorBarColor);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0.8);
    hold on;  
    plot(freq, GIRF_FT_mean_abs, 'b', 'LineWidth', 1);
    xlim(dispFreqRange); ylim([0, 1.1]);
    xlabel('Frequency [kHz]'); ylabel('Magnitude of GIRF [AU]');
    title('Magnitude of GIRF in Frequency Domain');
    hold off;
    
    subplot(212);
    plot(freq, GIRF_FT_mean_phase, 'b', 'LineWidth', 1);
    xlim(dispFreqRange); ylim([-4,4]);
    xlabel('Frequency [kHz]'); ylabel('Phase of GIRF [rad]');
    title('Phase of GIRF in Frequency Domain');
end

%% Plot GIRF in Time Domain
% Since we may truncated the GIRF in frequency domain,
% We need a new temporal axis
fullRORange = linspace(-max(roTime)/2, max(roTime)/2, size(GIRF_Temporal,1)); % in us

dispROTimeRange = [-50, 100]; % in unit of us

figure(777);
set(gcf,'color','white');
plot(fullRORange, real(GIRF_Temporal));
% xlim([min(fullRORange), max(fullRORange)]);
xlim(dispROTimeRange)
xlabel('Time [\mus]'); ylabel('Magnitude of GIRF [1/\mus]');
title('Magnitude of GIRF in Time Domain');
