% ************************************************************************
% Function: regDecomp
% Purpose:  Analyse quality of a registration
%
% Parameters:
%       XFd: original curves
%       XFdReg: registered curves
%       warpFd: warp curves
%
% Output:
%       decomp: registration decomposition structure
%
% ************************************************************************


function decomp = regDecomp( XFd, XFdReg, warpFd )

[ ampVar, phaVar, rSq, c ] = AmpPhaseDecomp( XFd, XFdReg, warpFd );

disp([ 'AmpVar = ' num2str( ampVar ) ...
       '; PhaVar = ' num2str( phaVar ) ...
       '; RSq = ' num2str( rSq ) ...
       '; C = ' num2str( c ) ]);

decomp.ampVar = ampVar;
decomp.phaVar = phaVar;
decomp.rSq = rSq;
decomp.c = c;


end