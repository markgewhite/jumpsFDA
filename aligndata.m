% ************************************************************************
% Function: aligndata
% Purpose:  Align time series data in preparation for FDA 
%           based on different methods
%
% Parameters:
%   x: cell array of time series
%   method: specified method either, 'Normalise','PadStart','PadEnd'
%   nReq: required length (for normalise)
%   tLength: number of points required
%   pad: padding value at start, if required
%
% Output:
%   ax: aligned time series as an array
%
%
% ************************************************************************


function ax = aligndata( x, method, nReq, tLength, pad, padMean )

if nargin < 6
    padMean = false;
end

nSubjects = size(x,1); % number of subjects
nTrials = size(x,2); % number of trials
nDim = size(x{1,1},2); % number of dimensions (assumes all are same)

switch method
    case 'Normalise'
        tspan = linspace(1,nReq,nReq); % define a normalised time series
        ax = zeros(nReq,nSubjects*nTrials,nDim); 
        for i = 1:nSubjects
            for j = 1:nTrials
                tLength = size(x{i,j},1);
                % normalise trial timespan to 0,1 for interpolation
                tLengthTrial = linspace(1,nReq,tLength);
                % interpolate to fit the standard timescale
                for k = 1:nDim
                    ax(:,(i-1)*nTrials+j,k) = ...
                            interp1(tLengthTrial,x{i,j}(:,k), ...
                                    tspan,'spline','extrap'); 
                end
            end
        end

    case 'PadStart'
        ax = zeros(tLength,nSubjects*nTrials,nDim);
        for i = 1:nSubjects
            for j = 1:nTrials
                tLengthTrial = min([size(x{i,j},1) tLength]);
                if padMean % use mean of initial values for padding values
                    pad = mean(x{i,j}(1:125,:),1);
                end
                for k = 1:nDim
                    ax(:,(i-1)*nTrials+j,k) = ...
                                [ ones(tLength-tLengthTrial,1)*pad(k); ...
                                    x{i,j}(end-tLengthTrial+1:end,k)];
                end
            end
        end
end

end

