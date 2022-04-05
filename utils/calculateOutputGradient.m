function [gradOutput, gradOutputFT] = calculateOutputGradient(sigS1, sigS2, refS1, refS2, params)
%Function for calculating the actual output gradient waveforms for GIRF calculation using the single polarity blips.
% function [gradOutput, gradOutputFT] = calculateOutputGradient(sigS1, sigS2, refS1, refS2, params)
% Algorithm: First perform temporal derivitive of the phase of raw T2* decay, then subtracte phases from two slices.
% 
% Input arguments:
%    sigS1 and sSigS2: raw T2* decays (with blips) from two symmetric slices
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

    nGradAmp = size(sigS1, 2); % Number of gradient blips
    nRep = size(sigS1, 3); % Number of repetition

    % Need to repeate if raw signal is not averaged
    % Note we are using mean of ref signal because they are independent to
    % the measurement of raw signal decay with blips; besides, the
    % repetition numbers are different in ref signal acquisition (5, vs 50 in raw signal)
    refS1 = repmat(refS1, [1,1,nRep]);
    refS2 = repmat(refS2, [1,1,nRep]);

    sigS1Corrected = zeros(size(sigS1), 'single');
    sigS2Corrected = zeros(size(sigS2), 'single');

    % Remember: We have ref scan before every three blip measurements
    % So nRefScan = floor((nGradAmp-1)/3)+1
    for nn = 1 : nGradAmp
        sigS1Corrected(:,nn,:) = sigS1(:,nn,:) ./ refS1(:,floor((nn-1)/3)+1,:);
        sigS2Corrected(:,nn,:) = sigS2(:,nn,:) ./ refS2(:,floor((nn-1)/3)+1,:);
    end

    sigS1Diff = diff(unwrap(angle(sigS1Corrected),[],1),1,1) / (params.adcDwellTime / 1e6); % in rad/s
    sigS2Diff = diff(unwrap(angle(sigS2Corrected),[],1),1,1) / (params.adcDwellTime / 1e6);

    % Compensate the one less point on readout due to diff function
    sigS1Diff = cat(1, zeros(1, size(sigS1Diff,2), size(sigS1Diff,3)), sigS1Diff);
    sigS2Diff = cat(1, zeros(1, size(sigS2Diff,2), size(sigS2Diff,3)), sigS2Diff);

    gradOutput = (sigS1Diff - sigS2Diff) / 2; % Phase averaged from two slices, in unit of rad/s
    gradOutput = gradOutput / (params.gammabar * 2*pi) / params.slicePos; % in unit of mT/m
    gradOutput = single(gradOutput);

    gradOutputFT = fftshift(fft(fftshift(gradOutput,1),[],1),1);

end