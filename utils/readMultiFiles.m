function varargout = readMultiFiles(varName, varargin)
%Function for loading variables with the same name from multiple files.
% function varargout = readMultiFiles(varName, varargin)
% 
% Input arguments:
%    The inputs are file names that are going to be loaded
% Output arguments:
%    Variables containing raw T2* decays from the corresponding files
% 
% Example:
%    To load three T2* decays from three files with variable name 'signal'
%    
%    [S1, S2, S3] = readMultiFiles('signal', 'FileName1', 'FileName2', FileName3');

% Author: Zhe "Tim" Wu
% Created: March 10, 2022

    if isempty(varargin)
        error("At least one file name/path is needed.");
    end
    
    if ~isa(varName, 'char')
        error("The first input arguments should be string for the variable name.");
    end
    
    if ~isa(varargin{1}, 'char')
        error("The second and later input arguments should be string for file paths.");
    end
    
    fn = varargin{1};
    s = load(fn);
    
    % Following parameters are the ones needs to be consistent between
    % files.
    roPts = s.roPts;
    nch = s.nch;
    dwellTime = s.dwellTime;
    xpts = size(s.(varName), 1);
    
    varargout{1} = s.(varName);
    clear s;
    
    for n = 2 : length(varargin)
        fn = varargin{n};
        s = load(fn);
        if s.roPts ~= roPts || s.nch ~= nch || s.dwellTime ~= dwellTime || size(s.(varName), 1) ~= xpts
            error("The numbers of readout points and coil channels, together with dwell time should be consistent between files.");
        end
        varargout{n} = s.(varName);
        clear s;
    end
end