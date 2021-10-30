% ************************************************************************
% Function: removeFaultyFd
% Purpose:  Remove extreme smoothed curves arising (it is assumed)
%           a faulty registration. They are identified by there being
%           an unreasonably large change in the area under the curve
%
% Parameters:
%       XFd: original smoothed curves
%       XFdReg: registered smoothed curves
%
% Output:
%       XFdTrim: smoothed curves without extremes
%       faulty: logical array identifying the curves removed
%
% ************************************************************************


function [ XFdTrim, XFdRegTrim, warpFdTrim, faulty ] = ...
                        removeFaultyFd( XFd, XFdReg, warpFd, setup )

% create a fine basis for integration
XBasis = getbasis( XFd );
XRng = getbasisrange( XBasis );
nXBasis = getnbasis( XBasis );

nFine = max( [201, 10*nXBasis] );
tFine = linspace( XRng(1), XRng(2), nFine )';

XPts = eval_fd( tFine, XFd );
XRegPts = eval_fd( tFine, XFdReg );
DWarpPts = eval_fd( tFine, warpFd, 1 );
D2WarpPts = eval_fd( tFine, warpFd, 2 );

% use area under the curve as criterion
XA = trapz( XPts );
XRegA = trapz( XRegPts );
DWarpA = trapz( DWarpPts );
D2WarpA = trapz( D2WarpPts );

AR = abs( XRegA./XA );
XZscore = abs(AR - mean(AR))/std(AR);

WZscore = abs(DWarpA-mean(DWarpA))/std(DWarpA);
W2Zscore = abs(D2WarpA-mean(D2WarpA))/std(D2WarpA);

figure(1);
plot( WZscore );
hold on;
plot( W2Zscore );
plot( XZscore );
hold off;

switch setup.faultCriterion
    case 'Warp'
        faulty = WZscore > setup.faultZScore ...
                            | W2Zscore > setup.faultZScore;
    case 'RelativeArea'
        faulty = XZscore > setup.faultZScore;
    otherwise
        error('Unrecognised fault criterion.');
end


% carry out removal using the coefficient matrix
coeff = getcoef( XFd );
coeff(:,faulty) = [];
XFdTrim = putcoef( XFd, coeff );

coeff = getcoef( XFdReg );
coeff(:,faulty) = [];
XFdRegTrim = putcoef( XFdReg, coeff );

coeff = getcoef( warpFd );
coeff(:,faulty) = [];
warpFdTrim = putcoef( warpFd, coeff );

disp(['Faulty registrations = ' num2str( sum(faulty) )]);

end
