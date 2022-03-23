% function varargout = readMultiFiles(varName, varargin)
% For loading raw GIRF data from multiple files
% Input arguments:
%    The inputs are file names that are going to be loaded
% Output arguments:
%    Variables containing raw T2* decays from the corresponding files
% 
% Created by Tim Wu, March 10, 2022

function varargout = readMultiFiles(varargin)
    if isempty(varargin)
        error("At least one file name/path is needed.");
    end
    
    if ~isa(varargin{1}, 'char')
        error("The second and later input arguments should be in characters for file paths.");
    end
    
    fn = varargin{1};
    s = load(fn);
    
    % Following parameters are the ones needs to be consistent between
    % files.
    roPts = s.roPts;
    nch = s.nch;
    dwellTime = s.dwellTime;
    xpts = size(s.kspace_all, 1);
    
    varargout{1} = s.kspace_all;
    clear s;
    
    for n = 2 : length(varargin)
        fn = varargin{n};
        s = load(fn);
        if s.roPts ~= roPts || s.nch ~= nch || s.dwellTime ~= dwellTime || size(s.kspace_all, 1) ~= xpts
            error("The numbers of readout points and coil channels, together with dwell time should be consistent between files.");
        end
        varargout{n} = s.kspace_all;
        clear s;
    end
end