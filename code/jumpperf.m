% ************************************************************************
% Function: jumpperf
% Purpose:  Calculate jump performance measures
%
% Parameters:
%       vgrf: discretised VGRF data
%
% Output:
%       perf.JHtov: jump height based on take-off velocity
%       perf.JHwd:  jump height based on work done
%       perf.PP:    peak power
%
% ************************************************************************


function perf = jumpperf( vgrf )

% initialise
g = 9.812; % acceleration due to gravity

N = size( vgrf, 2 );
perf.JHtov = zeros( N, 1);
perf.JHwd = zeros( N, 1);
perf.PP = zeros( N, 1);

for i = 1:N
    
    % calculate jump kinetics
    [ ~, vel, dis, pwr ] = kinetics( vgrf(:,i) );
   
    % jump height based on take-off velocity
    perf.JHtov(i) = 0.5*vel(end)^2/g;
    
    % jump height after adding in height at take-off
    perf.JHwd(i) = perf.JHtov(i) + dis(end);
    
    % peak power
    perf.PP(i) = max( pwr );
    
end

end


function [ a, v, s, p ] = kinetics( F )

g = 9.812; % acceleration due to gravity
rate = 1000; % sampling frequency

netF = F-1;

% calculate the acceleration
a = netF*g;

% calculate the velocity
v = cumtrapz( a )/rate;

% calculate the displacement
s = cumtrapz( v )/rate;

% calculate the power
p = F.*v*g;
       
end