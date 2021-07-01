% ************************************************************************
% Function: outputPCAtables
% Purpose:  Generate the PCA results tables (standard & long formats)
%           which will be saved to file
%
% Parameters:
%       measure: name of measure on which the analysis was done
%       ref: array of subject and jump number references in order of data
%       sData: array to translate subject number to subject ID
%       jumpOrder: lookup array for jump type
%       index: linking the subject to the jump order
%       pca: PCA functional data
%       oldperf: jump performance measure previously calculated
%       newperf: jump performance calculated with retained components
%
% Outputs:
%       results: table with results in standard format
%       results2: table with results in long format
%
% ************************************************************************


function [results,results2] = outputPCAtables( ...
                    measure,ref,sDataID,jumpOrder,index, ...
                    fd, pca,oldperf,newperf)

g = 9.812;
                
scoresFPC = pca.harmscr;
scoresACP = acp( fd, pca );
nJumps = size(scoresFPC,1);
nComponents = size(scoresFPC,2);
nVariables = size(scoresFPC,3);

nVar = 16;
numData = zeros(nVariables*nJumps,nVar+2*nComponents);
strData = strings(nVariables*nJumps,3);

% generate the output table

row = 0;
for v = 1:nVariables
    for w = 1:nJumps
        i = ref(w,1); % get the subject number (not ID)
        j = ref(w,2); % get the jump number (not ID)
        %row = (v-1)*nJumps+w;
        row = row+1;

        strData(row,1) = measure; % performance measure
        numData(row,1) = v; % variable counter
        numData(row,2) = sDataID(i); % subjectID
        numData(row,3) = j; % trial
        switch jumpOrder{sDataID(i)==index,j}
            case "V"
                strData(row,2) = "V"; % jumpType
                strData(row,3) = "WOA"; % arms condition
            case "VA"
                strData(row,2) = "V";
                strData(row,3) = "WA";
            case "H"
                strData(row,2) = "H";
                strData(row,3) = "WOA";
            case "HA"
                strData(row,2) = "H";
                strData(row,3) = "WA";
        end
        numData(row,4) = oldperf.height(i,j);
        numData(row,5) = newperf.height(w);
        numData(row,6) = oldperf.heightTO(i,j);
        numData(row,7) = oldperf.depth(i,j); 
        numData(row,8) = oldperf.peakPower(i,j)*g/oldperf.bwall(i);
        numData(row,9) = newperf.peakPower(w);
        numData(row,10) = oldperf.peakNegPower(i,j)*g/oldperf.bwall(i); 
        numData(row,11) = oldperf.workDone(i,j); 
        numData(row,12) = oldperf.time0(i,j); 
        numData(row,13) = oldperf.time1(i,j); 
        numData(row,14) = oldperf.time2(i,j); 
        numData(row,15) = oldperf.time3(i,j); 
        numData(row,16) = oldperf.bwall(i);
        for k = 1:nComponents
            numData(row,nVar+k) = round(scoresFPC(w,k,v),3); % PCA components
        end
        for k = 1:nComponents
            numData(row,nComponents+nVar+k) = round(scoresACP(w,k),3); % ACP components
        end
    end  
end
    
% create the output table by converting the data tables (strings in middle)

results = array2table(strData(:,1));
results.Properties.VariableNames{1} = 'Measure';

results = [results array2table(numData(:,1:3))];
results.Properties.VariableNames{2} = 'Variable';
results.Properties.VariableNames{3} = 'SubjectID';
results.Properties.VariableNames{4} = 'Trial';

results = [results array2table(strData(:,2:3))];
results.Properties.VariableNames{5} = 'JumpType';
results.Properties.VariableNames{6} = 'Arms';

results = [results array2table(numData(:,4:end))];
results.Properties.VariableNames{7} = 'Height';
results.Properties.VariableNames{8} = 'HeightRecalc';
results.Properties.VariableNames{9} = 'TakeOffHeight';
results.Properties.VariableNames{10} = 'CMD';
results.Properties.VariableNames{11} = 'PeakPower';
results.Properties.VariableNames{12} = 'PeakPowerRecalc';
results.Properties.VariableNames{13} = 'PeakNegPower';
results.Properties.VariableNames{14} = 'WorkDone';
results.Properties.VariableNames{15} = 'Time0';
results.Properties.VariableNames{16} = 'Time1';
results.Properties.VariableNames{17} = 'Time2';
results.Properties.VariableNames{18} = 'Time3';
results.Properties.VariableNames{19} = 'Bodyweight';
for i = 1:nComponents
    results.Properties.VariableNames{nVar+3+i} = char("PCA"+i);
end
for i = 1:nComponents
    results.Properties.VariableNames{nComponents+nVar+3+i} = char("ACP"+i);
end

% Create the long format table for SAS

for i = 1:nComponents
    batch = results(:,[1:nVar+2,nVar+2+i]);
    batch(:,20) = {i};
    batch = batch(:,[1:nVar+2,nVar+4,nVar+3]);
    batch.Properties.VariableNames{nVar+3} = 'PCA';
    % batch.Properties.VariableNames{nComponents+nVar+3} = 'ACP';
    batch.Properties.VariableNames{nVar+4} = 'Score';
    if i==1
        results2 = batch;
    else
        results2 = [results2; batch];
    end
end


end
