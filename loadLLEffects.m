% ************************************************************************
% Function: loadLLEffects
% Purpose:  Load parameters estimates from SAS's meta model (GLIMMIX)
%           a non-linear mixed model, too advanced for MATLAB to do
%           The file is a spreadsheet with three sheets.
%
% Parameters:
%
% Outputs:
%       data: table suitably structured for plotting
%
% ************************************************************************

function data = loadLLEffects( datapath, filename )

% read tables from the spreadsheet file
fullname = fullfile(datapath, filename);

importJH = readtable( fullname, ...
                detectImportOptions( fullname, 'Sheet', 'JH Params' ) );
importPP = readtable( fullname, ...
                detectImportOptions( fullname, 'Sheet', 'PP Params' ) );
importCL = readtable( fullname, ...
                detectImportOptions( fullname, 'Sheet', 'CL Params' ) );
            
% create the combined table
importJH.Properties.VariableNames{end} = 'JH_Estimate';
importPP.Properties.VariableNames{end} = 'PP_Estimate';
importCL.Properties.VariableNames{end} = 'CL_Estimate';

data = [ importJH importPP(:,end) importCL(:,end) ];


end

