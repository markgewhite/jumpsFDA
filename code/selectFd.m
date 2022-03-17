% ************************************************************************
% Function: selectFd
% Purpose:  Select a specified subset of functional curves
%
% Parameters:
%       Fd: smooth curves to select from
%       chosen: logical array of selected curves
%
% Output:
%       FdSelect: selected curves 
%
% ************************************************************************

function FdSelect = selectFd( Fd, chosen )

coeff = getcoef( Fd );
coeff( :, ~chosen ) = [];
FdSelect = putcoef( Fd, coeff );

end