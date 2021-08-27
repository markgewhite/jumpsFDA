% ************************************************************************
% Function: jumpInit
% Purpose:  Find jump initiation points
%
% Parameters:
%       X: VGRF data points
%       threshold1: primary jump detection threshold (large)
%       threshold2: secondary jump detection threshold (small)
%       sustainedPeriod: how long secondary threshold must be maintained
%
% Output:
%       t: jump initiation times (indices)
%
% ************************************************************************


function t = jumpInit( X, threshold1, threshold2, sustainedPeriod )

N = size( X, 2 );
t = zeros( N, 1 );

% jump initiation detection
for i = 1:N
    t0 = find( abs(X(:,i)-1) > threshold1, 1 );
    tx = find( abs(X(1:t0,i)-1) < threshold2, sustainedPeriod, 'last' )+1;
    t(i) = tx(1);
end

end
