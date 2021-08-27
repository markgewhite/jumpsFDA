% ************************************************************************
% Function: padData
% Purpose:  Pad time series to a specified length
%
% Parameters:
%       X: cell array of time series 
%       padLen: specified padded length
%       padValue: padding value
%
% Output:
%       Graphs
%
% ************************************************************************


function XP = padData( X, padLen, padValue )

% number of series
N = length( X );

% padded matrix of series
XP = zeros( padLen, N );

for i = 1:N
    % trial length
    trialLen = min( [ size(X{i}, 1), padLen] );
    % insert padding at beginning
    XP( :, i ) = [ ones( padLen-trialLen, 1 )*padValue; ...
                        X{i}(end - trialLen+1:end) ];
end
    
end