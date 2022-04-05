function [gResamp, roTime] = resampGradients(gradIn, params)
%Function for resampling input gradients into dwell time of ADC
% function [gResamp, roTime] = resampGradients(gradIn, params)
% 
% Input arguments:
%    gradIn: input gradients, in size of gradInNumber × gradAmplitudeNumbers
%    params.gradRasterTime: scalar, in unit of us
%    params.adcDwellTime: scalar, in unit of us
%    params.roPts: scalar, number of adc samples
%    
% Output arguments:
%    gResamp: resampled gradients, in size of adcSampNum × gradAmplitudeNumbers
%    roTime: time for each resampled gradient, in size of adcSampNum × 1
% 
% Author: Zhe "Tim" Wu
% Created: March 10, 2022

    nGradSamples = size(gradIn, 1);
    gradTime = 0:(nGradSamples - 1);
    gradTime = gradTime * params.gradRasterTime;
    gradTime = gradTime(:); % Make it as a column vector

    roTime = 0:(params.roPts - 1);
    roTime = roTime * params.adcDwellTime;
    roTime = roTime(:); % Make it as a column vector

    nGradIn = size(gradIn,2); % Number of input gradients

    % We have already have the timeline for raw signal as roTime
    gResamp = zeros(params.roPts, nGradIn, 'single');
    for n = 1 : nGradIn
        gResamp(:, n) = interp1(gradTime, gradIn(:,n), roTime);
    end
    gResamp(end,:) = 0; % Fix the possible interp error on the last point.

end