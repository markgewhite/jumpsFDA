% ************************************************************************
% Function: encodeProc
% Purpose:  Encode the processing procedure
%
% Parameters:
%           i, j, k, l : loop counters
%
% Output:
%           name : string encoding
%           regCode: : registration encoding
%
% ************************************************************************

function [ name, regCode ] = encodeProc( i, j, k, l )

kbin = dec2bin( k-1, 4 );

if l==1
    cbin = '-';
else
    cbin = 'C';
end

regCode.grfmin = (kbin(1)=='1');
regCode.pwrmin = (kbin(2)=='1');
regCode.pwrcross = (kbin(3)=='1');
regCode.pwrmax = (kbin(4)=='1');
regCode.ct = (l==2);

name = [ num2str(i) '-' num2str(j) 'VGRF' kbin cbin ];

end