% ************************************************************************
% Function: findGRFlandmarks
% Purpose:  Locate landmarks for a GRF set of curves
%
% Parameters:
%       tspan: timespan for the GRF series
%       grf: GRF functional data object
%       lm: landmarks structure including:
%          .grfmin   % use VGRF minimum as a landmark?
%          .pwrmin   % use Power minimum as a landmark?
%          .pwrcross % use Power crossing point as a landmark?
%          .pwrmax   % use Power maximum as a landmark?e-off
%
% Output:
%       landmarks: where the landmarks are located in time (means & cases)
%
% ************************************************************************


function landmarks = findGRFlandmarks( tspan, grfFd, lm )


% evaluate the points for the GRF curves and derivatives
% take off so the crossing point really is zero
grf = eval_fd( tspan, grfFd ) - 1;

% the discrete representation must be well behaved so 
% central differences can be used (no need for the denominator)
grfD1 = grf( 3:end, : ) - grf( 1:end-2, : );

% calculate the power curves (add back 1)
vel = cumtrapz( grf );
pwr = (grf+1).*vel;
% evaluate the points for the power curves and derivatives
pwrD1 = pwr( 3:end, : ) - pwr( 1:end-2, : );

n = size( grf, 2 ); % number of jumps

% initialise
tGrfMin = zeros( n, 1 );
tPwrMin = zeros( n, 1 );
tPwrCross = zeros( n, 1 );
tPwrMax = zeros( n, 1 );

for i = 1:n
  
    % --- Find GRF Minimum ---
    % find the time indices where there are minima, below zero
    % exclude the last 3% before take-off
    tend = floor( 0.99*length(grf) );
    grfD1Up = grfD1( 2:tend, i );
    grfD1Dn = grfD1( 1:tend-1, i );
    index = find( grfD1Dn.*grfD1Up<0 & grfD1Up>0 & grf(2:tend,i)<0 )+1;
    if isempty( index )
        index = 1;
    end
    % select the global minimum
    [ ~, minIndex0 ] = min( grf( index, i ) );
    grfMinIdx = index( minIndex0 );
    tGrfMin( i ) = tspan( grfMinIdx );
    
    % --- Find PWR Minimum ---
    % find the time indices where there are minima, below zero
    pwrD1Up = pwrD1( 2:tend, i );
    pwrD1Dn = pwrD1( 1:tend-1, i );
    index = find( pwrD1Dn.*pwrD1Up<0 & pwrD1Up>0 & pwr(2:tend,i)<0 )+1;
    if isempty( index )
        index = 1;
    end
    % select the global minimum
    [ ~, minIndex0 ] = min( pwr( index, i) );
    pwrMinIdx = index( minIndex0 );
    tPwrMin( i ) = tspan( pwrMinIdx );
    
    % --- Find PWR Crossover ---
    % find the t-axis crossing indices in PWR derivative after the minimum
    pwrD0Up = pwr( pwrMinIdx+1:tend, i );
    pwrD0Dn = pwr( pwrMinIdx:tend-1, i );   
    index = find( pwrD0Dn.*pwrD0Up<0 & pwrD0Up>0 ) ;
    if isempty( index )
        index = 10;
    end
    % select the last crossing point as the minimum GRF
    pwrCrossIdx = index(end) + pwrMinIdx;
    tPwrCross( i ) = tspan( pwrCrossIdx );
    
    % --- Find PWR Maximum ---
    % find the t-axis crossing indices in PWR derivative after crossover
    index = find( pwrD1Dn.*grfD1Up<0 & pwrD1Up<0 & pwr(2:tend,i)>0 )+1;
    % select the last crossing point as the minimum GRF
    if isempty( index )
        index = tend;
    end
    % select the global maximum
    [ ~, maxIndex0 ] = max( pwr( index, i ) );
    pwrMaxIdx = index( maxIndex0 );
    tPwrMax( i ) = tspan( pwrMaxIdx );
    
  
end

% take means
tGrfMinMean = mean( tGrfMin );
tPwrMinMean = mean( tPwrMin );
tPwrCrossMean = mean( tPwrCross );
tPwrMaxMean = mean( tPwrMax );

% assemble return array structure

landmarks.mean = 0;
landmarks.case = zeros( n, 1 );

if lm.grfmin
    landmarks.mean = [ landmarks.mean, tGrfMinMean ];
    landmarks.case = [ landmarks.case, tGrfMin ];
end
if lm.pwrmin
    landmarks.mean = [ landmarks.mean, tPwrMinMean ];
    landmarks.case = [ landmarks.case, tPwrMin ];
end
if lm.pwrcross
    landmarks.mean = [ landmarks.mean, tPwrCrossMean ];
    landmarks.case = [ landmarks.case, tPwrCross ];
end
if lm.pwrmax
    landmarks.mean = [ landmarks.mean, tPwrMaxMean ];
    landmarks.case = [ landmarks.case, tPwrMax ];
end

landmarks.mean(1) = [];
landmarks.case(:,1) = [];


end


