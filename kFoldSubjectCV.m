% ************************************************************************
% Function: kFoldSubjectCV
% Purpose:  Define a set of K-Fold partitions by subject
%
% Parameters:
%       subjectID: list of subject identifiers by case
%       nRepeats: number of CV repeats
%       nFolds: number of folds
%
% Outputs:
%       trnSelect: training partitions
%       tstSelect: testing partitions (inverse)
%
% ************************************************************************

function [trnSelect, tstSelect] = kFoldSubjectCV( subjects, nRepeats, nFolds ) 

% Repeated K-fold cross validation by subject

trnSelect = false( length( subjects ), nRepeats*nFolds );
tstSelect = false( length( subjects ), nRepeats*nFolds );

% identify unique subjects
sID = unique( subjects );

% for reproducibility across data sets

for i = 1:nRepeats
    
    % partition based on subjects
    subjectCV = cvpartition( length(sID), 'KFold', nFolds );

    for k = 1:nFolds
        % extract the logical array of testing and training
        trnSelect( :, (i-1)*nFolds+k ) = ismember( subjects, ...
                                sID(subjectCV.training(k)) );
        tstSelect( :, (i-1)*nFolds+k ) = ismember( subjects, ...
                                sID(subjectCV.test(k)) );
    end
        
end

end