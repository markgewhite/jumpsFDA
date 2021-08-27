% ************************************************************************
% Function: timeNormVGRF
% Purpose:  Time normalise the VGRF, starting at jump initiation
%
% Parameters:
%       XFd: smoothed VGRF curves
%       setup: relevant settings, including jump detection threshold
%
% Output:
%       XNFd: time-normalised smoothed VGRF curves
%
% ************************************************************************


function XNFd = timeNormVGRF( XFd, t, setup )

% convert smooth curve to time series 
X = eval_fd( t, XFd );

N = size( X, 2 );
initTime = zeros( N, 1 );

% jump initiation detection
for i = 1:N
    initTime(i) = find( abs(X(:,i)-1) > setup.detectThreshold, 1 );   
end

% median length gives standard length
fixLen = size( X, 1 ) - median( initTime );

% time normalise to this length
XN = aligndata( X, 'Normalise', ...
                    setup.tFreq, ...
                    fixLen, ...
                    setup.initial );

% set new time span
tn = -fixLen+1:0;

% smooth the time normalised data
XNFd = smoothVGRF( XN, tn, setup );


end
