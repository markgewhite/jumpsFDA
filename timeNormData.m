% ************************************************************************
% Function: timeNormData
% Purpose:  Time normalise data to list of specified lengths.
%           The time series have already been padded to a fixed length.
%
% Parameters:
%       X: cell array of time series
%       fixedLen: fixed length for time normalisation
%
% Output:
%       XN: time normalised time series
%
% ************************************************************************


function XN = timeNormData( X, fixedLen )

% number of series
N = length( X );

% padded matrix of series
XN = zeros( fixedLen, N );

% define standard time series
tspan1 = linspace( -fixedLen+1, 0, fixedLen ); 

for i = 1:N
    
    % define current time series
    tspan0 = linspace( -fixedLen+1, 0, length(X{i}) );
    
    % interpolate to fit the standard timescale
    XN( :, i ) = interp1( tspan0, ...
                          X{i}, ...
                          tspan1, ...
                          'spline', 'extrap' ); 
end

    
end