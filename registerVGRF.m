% ************************************************************************
% Function: registerVGRF
% Purpose:  Perform curve registration on VGRF data
%
% Parameters:
%       t: time span
%       XFd: smoothed curves
%       setup: smoothing settings
%       opt: relevant options
%
% Output:
%       XFdReg: registered smoothed curves
%
% ************************************************************************


function [ XFdReg, warpFd, decomp ] = registerVGRF( t, XFd, setup, opt )

% initialise
monotonic = true;
nProcustes = setup.nIterations;
N = size( getcoef( XFd ), 2 );

XFdReg = XFd;
warpT = repmat( t, N, 1 )';

% plot curves
%subplot( N+1, 2, 1 );
%plot( XFdReg );
%drawnow;

% use a Procustes style loop
for i = 1:nProcustes
    
    % locate landmarks
    lm = findGRFlandmarks( t, XFdReg, setup.lm, opt );
    
    % sort into order (in necessary)
    [ lm.mean, lm.order ] = sort( lm.mean, 'ascend' );
    lm.case = lm.case( :, lm.order );
    
    if isempty( lm.mean )
        
        % *** continuous registration  ***
        
        wBasis = create_bspline_basis( [t(1),t(end)], ...
                                       setup.nBasis, ...
                                       setup.basisOrder );
        
        wFdReg = fd( zeros( setup.nBasis, 1), wBasis );
        wFdRegPar = fdPar( wFdReg, 1, setup.wLambda );
                              
        XMeanFd = mean( XFdReg );
        
        [ XFdReg, warpFd ] = register_fd( ...
                                            XMeanFd, XFdReg, wFdRegPar) ;                       
                                   
    else
        
        % *** landmark registration ***
        
        disp(['Landmark means = ' num2str( lm.mean )]);
        disp(['Landmark SDs   = ' num2str( std( lm.case ) )]);
                        
        wBasis = create_bspline_basis( [t(1),t(end)], ...
                                       setup.nBasis, ...
                                       setup.basisOrder, ...
                                       [ t(1) lm.mean t(end) ] );

        wFdReg = fd( zeros( setup.nBasis, 1), wBasis );
        wFdRegPar = fdPar( wFdReg, 1, setup.wLambda );

        [ XFdReg, warpFd ] = landmarkreg( ...
                                    XFdReg, lm.case, lm.mean, ...
                                    wFdRegPar, monotonic, setup.XLambda );
                            
    end
    
    % update time warp to maintain link with original
    for j = 1:N
        % warp the previous warp
        warpT(:,j) = eval_fd( warpT(:,j), warpFd(j) );
        % separate the points evenly
        warpT(:,j) = interp1( t, warpT(:,j), t, 'spline', 'extrap' );
        % impose limits in case of over/underflow
        warpT(:,j) = max( min( warpT(:,j), t(end) ), t(1) );
    end
    warpFd = smooth_basis( t, warpT, wFdRegPar );
    
    % compute amplitude and phase decomposition
    [ ampVar, phaVar, rSq, c ] = AmpPhaseDecomp( XFd, XFdReg, warpFd );

    disp([ 'Iteration = ' num2str(i) ...
           '; AmpVar = ' num2str( ampVar ) ...
           '; PhaVar = ' num2str( phaVar ) ...
           '; RSq = ' num2str( rSq ) ...
           '; C = ' num2str( c ) ]);
    
    %subplot( N+1, 2, i*2+1 );
    %plot( XFdReg );
    %subplot( N+1, 2, i*2+2 );
    %plot( warpFd );
    %drawnow;
    
end

if ~isempty( lm.mean )
    % locate landmarks to determine end point
    lm = findGRFlandmarks( t, XFdReg, setup.lm, opt );
    disp(['Final Landmark means = ' num2str( lm.mean )]);
    disp(['Final Landmark SDs   = ' num2str( std( lm.case ) )]);
end
    
    
decomp.ampVar = ampVar;
decomp.phaVar = phaVar;
decomp.rSq = rSq;
decomp.c = c;


end