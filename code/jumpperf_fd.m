% ************************************************************************
% Function: jumpperf_fd
% Purpose:  Calculate jump performance measures from smoothed curves
%
% Parameters:
%       vgrfFd: smoothed VGRF curves
%
% Output:
%       perf.JHtov: jump height based on take-off velocity
%       perf.JHwd:  jump height based on work done
%       perf.PP:    peak power
%
% ************************************************************************


function perf = jumpperf_fd( vgrfFd )

% set time span for evaluation
tRng = getbasisrange( getbasis(vgrfFd) );
tSpan = tRng(1):tRng(2);

% discretise all curves
vgrf = eval_fd( tSpan, vgrfFd ); 

perf = jumpperf( vgrf );

end
