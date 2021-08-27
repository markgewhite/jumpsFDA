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


function ax = aligndata( X, method, nReq, tLength, pad, padMean )

if nargin < 6
    padMean = false;
end

nSubjects = size( X, 1 ); % number of subjects
nTrials = size( X, 2 ); % number of trials

switch method
    
    case 'Normalise'
        tspan = linspace( 1, nReq, nReq ); % define a normalised time series
        ax = zeros( nReq, nSubjects*nTrials ); 
        for i = 1:N
            tLength = length( X(:, 1 );
                % normalise trial timespan to 0,1 for interpolation
                tLengthTrial = linspace( 1, nReq, tLength );
                % interpolate to fit the standard timescale
                ax( :, (i-1)*nTrials+j ) = ...
                            interp1(tLengthTrial,X{i,j}(:,, ...
                                    tspan,'spline','extrap'); 
        end

    case 'PadStart'
        ax = zeros(tLength,nSubjects*nTrials,nDim);
        for i = 1:nSubjects
            for j = 1:nTrials
                tLengthTrial = min([size(X{i,j},1) tLength]);
                if padMean % use mean of initial values for padding values
                    pad = mean(X{i,j}(1:125,:),1);
                end
                for k = 1:nDim
                    ax(:,(i-1)*nTrials+j,k) = ...
                                [ ones(tLength-tLengthTrial,1)*pad(k); ...
                                    X{i,j}(end-tLengthTrial+1:end,k)];
                end
            end
        end
end

end

