% ************************************************************************
% Function: timeNormData
% Purpose:  Time normalise data to list of specified lengths.
%           The time series have already been padded to a fixed length.
%
% Parameters:
%       X: matrix of time series
%       trialLen: specified trial lengths
%       fixedLen: fixed length for time normalisation
%
% Output:
%       XN: time normalised time series
%
% ************************************************************************


function XN = timeNormData( X, trialStarts, fixStart )

% number of series
N = size( X, 2 );

% time normalised length
fixLen = size( X, 1 ) - fixStart + 1;

% trial lengths
trialLen = size( X, 1 ) - trialStarts + 1;

% padded matrix of series
XN = zeros( fixLen, N );

% define standard time series
tspan1 = linspace( -fixLen+1, 0, fixLen ); 

for i = 1:N
    
    % define current time series
    tspan0 = linspace( -fixLen+1, 0, trialLen(i) );
    
    % interpolate to fit the standard timescale
    XN( :, i ) = interp1( tspan0, ...
                          X( trialStarts(i):end, i ), ...
                          tspan1, ...
                          'spline', 'extrap' ); 
end

    
end