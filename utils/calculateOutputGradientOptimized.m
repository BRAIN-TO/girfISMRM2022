function [gradOutput, gradOutputFT] = calculateOutputGradientOptimized(sigS1_POS, sigS1_NEG, sigS2_POS, sigS2_NEG, refS1, refS2, params)
%Function calculating the actual output gradient waveforms for GIRF calculation using the dual polarity (optimized) method.
% function [gradOutput, gradOutputFT] = calculateOutputGradientOptimized(sigS1_POS, sigS1_NEG, sigS2_POS, sigS2_NEG, refS1, refS2, params)
% First perform temporal derivitive of the phase of raw T2* decay, then subtracte phases from two slices.
% 
% Input arguments:
%    sigS1 and sSigS2: raw T2* decays (with blips) from two symmetric slices
%                           The suffix POS/NEG means the blip is applied in
%                           positive or negative polarities.
%                           in size of roPts × nGradAmplitudes × nRepetition
%    refS1 and refS2: raw T2* decays (without blips) from two symmetric slices
%                           in size of roPts × nGradAmplitudes
%    params: parameters for calculations. A structure.
%    
% Output arguments:
%    gradOutput: actual output gradients in time domain with the same size as inputs
%    gradOutputFT: actual output gradients in frequency domain with the same size as inputs
% 
% Author: Zhe "Tim" Wu
% Created: March 10, 2022

    nGradAmp = size(sigS1_POS, 2); % Number of gradient blips
    nRep = size(sigS1_POS, 3); % Number of repetition

    % Need to repeate if raw signal is not averaged
    % Note we are using mean of ref signal because they are independent to
    % the measurement of raw signal decay with blips; besides, the
    % repetition numbers are different in ref signal acquisition (5, vs 50 in raw signal)
    refS1 = repmat(refS1, [1,1,nRep]);
    refS2 = repmat(refS2, [1,1,nRep]);

    sigS1Corrected_POS = zeros(size(sigS1_POS), 'single');
    sigS1Corrected_NEG = zeros(size(sigS1_NEG), 'single');
    sigS2Corrected_POS = zeros(size(sigS2_POS), 'single');
    sigS2Corrected_NEG = zeros(size(sigS2_NEG), 'single');

    % Remember: We have ref scan before every three blip measurements
    % So nRefScan = floor((nGradAmp-1)/3)+1
    for nn = 1 : nGradAmp
        sigS1Corrected_POS(:,nn,:) = sigS1_POS(:,nn,:) ./ refS1(:,floor((nn-1)/3)+1,:);
        sigS1Corrected_NEG(:,nn,:) = sigS1_NEG(:,nn,:) ./ refS1(:,floor((nn-1)/3)+1,:);
        sigS2Corrected_POS(:,nn,:) = sigS2_POS(:,nn,:) ./ refS2(:,floor((nn-1)/3)+1,:);
        sigS2Corrected_NEG(:,nn,:) = sigS2_NEG(:,nn,:) ./ refS2(:,floor((nn-1)/3)+1,:);
    end

    sigS1Diff_POS = diff(unwrap(angle(sigS1Corrected_POS),[],1),1,1) / (params.adcDwellTime / 1e6); % in rad/s
    sigS1Diff_NEG = diff(unwrap(angle(sigS1Corrected_NEG),[],1),1,1) / (params.adcDwellTime / 1e6);
    sigS2Diff_POS = diff(unwrap(angle(sigS2Corrected_POS),[],1),1,1) / (params.adcDwellTime / 1e6);
    sigS2Diff_NEG = diff(unwrap(angle(sigS2Corrected_NEG),[],1),1,1) / (params.adcDwellTime / 1e6);

    % Compensate the one less point on readout due to diff function
    sigS1Diff_POS = cat(1, zeros(1, size(sigS1Diff_POS,2), size(sigS1Diff_POS,3)), sigS1Diff_POS);
    sigS1Diff_NEG = cat(1, zeros(1, size(sigS1Diff_NEG,2), size(sigS1Diff_NEG,3)), sigS1Diff_NEG);
    sigS2Diff_POS = cat(1, zeros(1, size(sigS2Diff_POS,2), size(sigS2Diff_POS,3)), sigS2Diff_POS);
    sigS2Diff_NEG = cat(1, zeros(1, size(sigS2Diff_NEG,2), size(sigS2Diff_NEG,3)), sigS2Diff_NEG);

    gradOutput = ((sigS1Diff_POS - sigS1Diff_NEG) + (sigS2Diff_NEG - sigS2Diff_POS)) / 4; % Phase averaged from two slices with both positive/negative blips, in unit of rad/s;
    gradOutput = gradOutput / (params.gammabar * 2*pi) / params.slicePos; % in unit of mT/m
    gradOutput = single(gradOutput);

    gradOutputFT = fftshift(fft(fftshift(gradOutput,1),[],1),1);

end