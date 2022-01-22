% ************************************************************************
% Function: jumpFlightTime
% Purpose:  Find the the jump flight time based on the vGRF 
%           Flight is indicated by where (vGRF < 10 N)
%
% Parameters:
%   vgrf: array with vertical GRF
%   rate: sampling rate used
%
% Output:
%    ftime: flight time
%
% ************************************************************************


function ftime = jumpFlightTime( VGRF, rate )

% constants
minGRF = 10; % Newtons

% find first point where vGRF falls below minimum threshold (take-off)
t1 = find(abs(VGRF)<minGRF,1);
if isempty(t1)
    t1 = 1;
end

% find first point after t1 where vGRF rises above min threshold (landing)
% allow for oscillation of the force platform with another 10 time intervals 
t2 = find(abs(VGRF(t1+10:end))>minGRF,1)+t1;
if isempty(t2)
    t2 = size(VGRF,1);
end

% calculate flight time
ftime = (t2-t1)/rate;

end
