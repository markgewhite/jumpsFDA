% ************************************************************************
% Function: performFDA
% Purpose:  Carry out a series of functional data analyses on a dataset
%           that includes time normalisation, filtering, smoothing, 
%           registration and principal components analysis
%
% Parameters:
%       x: the dataset in question
%       timespan: time domain
%       takeoff: takeoff (index)
%       setup: list of settings to use
%       option: list of processing options
%       LMfunc: function handle for locating landmarks
%       selectID: selection of curves for mean (optional)
%
% Output:
%       xFd: smoothed functional data object
%       xFdPar: functional data object's parameters
%       pca: sets of principal components analysis
%       warpFd: timewarped functional data object
%       landmarks: where the landmarks are located in time (means & cases)
%       xDecomp: results vector of amplitude/phase variance decomposition
%
% ************************************************************************


function [ xFd, xFdPar, xPCA, warpFd, landmarks, xDecomp ] = ...
                        performFDA( x, timespan, takeoff, ...
                                    setup, option, ...
                                    LMfunc, meanID )


% ************************************************************************
%   Constants
% ************************************************************************

nJumps = size(x,1);
tFreq = setup.data.tFreq;
maxLength = setup.data.maxLength;

warpFd = [];
landmarks = [];
xDecomp = zeros(1,6);


% ************************************************************************
%   Normalise time
% ************************************************************************

if option.doTimeNormalisation
    xNorm = aligndata( x, 'Normalise', maxLength );
else
    if option.doTruncation
        xNorm = aligndata(  x, 'PadStart', ...
                            setup.data.tFreq, ...
                            setup.data.maxLength, ...
                            setup.data.initial );
    else
        xNorm = pad2phase(  x, takeoff, ...
                            setup.data.maxLength1, ...
                            setup.data.maxLength2 );
    end
end
tStart = timespan(1);
tEnd = timespan(end);


% ************************************************************************
%   Filter data
% ************************************************************************

if option.doFiltering
    filterSpec = designfilt('lowpassfir','FilterOrder',4, ...
                        'CutoffFrequency',setup.data.cutoffFreq,'SampleRate',tFreq);
    [butterB,butterA] = butter(4,setup.data.cutoffFreq/(tFreq/2),'low');
    xPadRange = 500; % buffering
    xPadded = [ones(xPadRange,nJumps); xNorm; zeros(xPadRange,nJumps)];
    xFilt = zeros(size(xPadded));
    xFiltB = zeros(size(xPadded));
    sse = zeros(nJumps,1);
    % delay = mean(grpdelay(filterSpec,tNorm+2*xPadRange,tNorm+2*xPadRange));
    lag = 28; % phase lag
    for i = 1:nJumps
        xFilt(:,i) = filter(filterSpec,xPadded(:,i));
        xFilt(:,i) = [xFilt(lag+1:end,i); zeros(lag,1)];
        xFiltB(:,i) = filter(butterB,butterA,xPadded(:,i));
        xFiltB(:,i) = [xFiltB(lag+1:end,i); zeros(lag,1)];
        sse(i) = sum((xFiltB(xPadRange+1:xPadRange+maxLength,i)-xNorm(:,i)).^2);
    end

    stderr = sqrt(sum(sse))/(nJumps*(maxLength-4));
    disp(['Standard error from low-pass filter = ',num2str(stderr) ])

    xNorm = xFilt(xPadRange+1:xPadRange+maxLength,:);
    
    disp(['RMSE from filtering = ' ...
                    num2str(perf_rmse(xNorm,xCrit,tFreq,1)) ' m']);
    
end


% ************************************************************************
%   Setup basis
% ************************************************************************

% set up a basis and fdPar object for smoothing the data
basis = create_bspline_basis([tStart,tEnd],setup.Fd.basis,setup.Fd.basisOrder);
xFdPar = fdPar(basis,setup.Fd.penaltyOrder,setup.Fd.lambda); 





% ************************************************************************
%   Create FD object
% ************************************************************************

[xFd,df,gcv,~,ssErr] = smooth_basis(timespan,xNorm,xFdPar);
xFd = putnames(xFd,setup.Fd.names);

disp(['Degrees of freedom = ',num2str(df)])
disp(['Generalized cross-validation = ',num2str(sum(gcv,2)')])

%  estimate standard error of fit
stderr = sqrt(sum(ssErr,2)/(nJumps*(maxLength-df))); 
disp(['Standard error of fit = ',num2str(stderr)])


% ************************************************************************
%   Check the goodness of fit
% ************************************************************************

if option.doCheckFit
    %  plot each curve along with the data
    clf;
    plotfit_fd(xNorm,timespan,xFd,option.doResidualCheck,option.doSortCheck);
end


% ************************************************************************
%   Locate the landmarks
% ************************************************************************

if option.lm.doReg0 || option.lm.doReg1
    % produce a heavily smoothed curve for finding landmarks only
    xFdParLMSearch = fdPar( basis, setup.Fd.penaltyOrder, ...
                                 setup.lm.searchLambda );
    xFdLMSearch = smooth_basis( timespan, xNorm, xFdParLMSearch );
    
    landmarks = LMfunc( timespan, xFdLMSearch, setup.lm.set, option.lm );
    outliers = find( abs(landmarks.case-landmarks.mean)>std(landmarks.case)*3 );
    disp(['Number Outlier Landmarks = ' num2str(length(outliers))]);
end


% ************************************************************************
%   Carry out landmarks registration (non-derivative)
% ************************************************************************

if option.lm.doReg0
    %  carry out the landmark registration
    if option.lm.doReg0Calc
        % do the intensive calculations
     
        % add an extra basis for extra landmarks
        nBasisLM = setup.lm.basis+length(landmarks.mean)-1;
        
        [ xFdLM, warpFd, wFd ] = curveregistration( ...
                                timespan,xFd,[],landmarks, ...
                                nBasisLM, setup.lm.basisOrder, ...
                                setup.lm.lambda, setup.lm.yLambda);  
        save(setup.lm.filename,'xFdLM','warpFd','wFd');
    else
        % read from a file
        load(setup.lm.filename);
    end
    
    % Find out how much change has been introduced by registration
    % as indicated by the correlation between before and after
    [ ampVar, phaVar, rSq, c ] = ...
                AmpPhaseDecomp( xFd, xFdLM, warpFd );
    xDecomp(1) = ampVar;
    xDecomp(2) = phaVar;
    xDecomp(3) = rSq;
    xDecomp(4) = c;
    
    rmse0 = mean(sqrt(sum(eval_fd(timespan, xFd-mean(xFd)).^2,2))/length(timespan));
    rmse1 = mean(sqrt(sum(eval_fd(timespan, xFdLM-mean(xFdLM)).^2,2))/length(timespan));
    xDecomp(5) = rmse1;
    xDecomp(6) = rmse0;
    
    disp(['Decomposition: ' num2str(xDecomp)]);
    xFdLMSearchPts = eval_fd( timespan, xFdLM );
    xFdLMSearch = smooth_basis( timespan, xFdLMSearchPts, xFdParLMSearch );
    landmarks2 = LMfunc( timespan, xFdLMSearch, setup.lm.set, option.lm );
    disp(['Landmark SD: ' num2str( std( landmarks2.case, 1 ) ) ]);
    
    xFd = xFdLM;
    
    %plot(xFdLM);
    %drawnow;
    %plot(warpFd);
    %pause;
    %plot(wFd);
    %pause;

end


% ************************************************************************
%   Perform continuous registration
% ************************************************************************

if option.ct.doReg
    
    % carry out the registration
    if option.ct.doRegCalc
        % do the intensive calculations
        [ xFdCT, warpFd, wFd ] = curveregistration( ...
                                        timespan, xFd, [], [], ...
                                        setup.ct.basis, ...
                                        setup.ct.basisOrder, ...
                                        setup.ct.lambda);  
        save(setup.ct.filename,'xFdCT','warpFd','wFd');
    else
        load(setup.ct.filename);
    end
    
    % Find out how much change has been introduced by registration
    % as indicated by the correlation between before and after
    [ ampVar, phaVar, rSq, c ] = ...
                AmpPhaseDecomp( xFd, xFdCT, warpFd );
    xDecomp(1) = ampVar;
    xDecomp(2) = phaVar;
    xDecomp(3) = rSq;
    xDecomp(4) = c;
        
    rmse0 = mean(sqrt(sum(eval_fd(timespan, xFd-mean(xFd)).^2,2))/length(timespan));
    rmse1 = mean(sqrt(sum(eval_fd(timespan, xFdCT-mean(xFdCT)).^2,2))/length(timespan));
    xDecomp(5) = rmse1;
    xDecomp(6) = rmse0;
    
    xFd = xFdCT;
    
    plot(xFdCT);
    drawnow;
    
end

                
% ************************************************************************
%   Perform principal components analysis
% ************************************************************************

if ~isempty(meanID)
    % use a different mean for the reference curve
    xrefFd = mean(xFd(meanID));
    xPCA.unrotated = pca_fd_setmean( xFd, setup.pca.nComponents, ...
                                                    xFdPar, 1, xrefFd );
else
    % use the overall mean for the reference curve
    xPCA.unrotated = pca_fd(    xFd, ...
                                setup.pca.nComponents, ...
                                xFdPar, ...
                                option.pca.doCentreFunctions );
end

if option.lm.doReg0 || option.lm.doReg1 || option.ct.doReg
    xPCA.warp = pca_fd(     warpFd, ...
                            setup.pca.nComponents, ...
                            xFdPar, ...
                            option.pca.doCentreFunctions );
    xPCA.time = pca_fd(     wFd, ...
                            setup.pca.nComponents, ...
                            xFdPar, ...
                            option.pca.doCentreFunctions );
end

disp(['Total proportion of variation accounted for = ' ...
                        num2str(sum(xPCA.unrotated.varprop))]);
   
if option.pca.doVarimaxRotation
    xPCA.varimax = varmx_pca(xPCA.unrotated);
    if option.lm.doReg0 || option.lm.doReg1 || option.ct.doReg
        xPCA.warpVarimax = varmx_pca(xPCA.warp);
        xPCA.timeVarimax = varmx_pca(xPCA.time);
    end
end

if option.pca.doShowComponents
    clf;
    plot_pca_fd(xPCA.unrotated);
    if option.pca.doVarimaxRotation
        plot_pca_fd(xPCA.varimax);
    end
end
    

% ************************************************************************
%   Validate how many components to retain based on jump performance
% ************************************************************************

if option.pca.doComponentsValidation
    
    rmsei = zeros(setup.pca.nComponents,1);
    is1D = (size(xNorm,3)==1);
    drop = -0.25;
    for i = 1:setup.pca.nComponents
        if is1D
            pca_perf = jumpheight_pca(timespan,xPCA.unrotated,i);
        else
            pca_perf = jumpdistance_pca(timespan,xPCA.unrotated,[],drop,i);
        end
        crit_perf = zeros(nJumps,1);
        if is1D
            for j = 1:nJumps
                crit_perf(j) = jumpheight(xNorm(:,j),tFreq,1);
            end
        else
            xPts = zeros(size(xNorm,1),2);
            for j = 1:nJumps
                % convert to Cartesian
                [xPts(:,1),xPts(:,2)] = pol2cart(xNorm(:,j,2),xNorm(:,j,1)); 
                % calculate jump distances
                crit_perf(j) = jumpdist(xPts,tFreq,1,drop);
            end
        end
        rmsei(i) = sqrt(sum((pca_perf-crit_perf).^2)/nJumps);
    end
    %  plot the results
    figure(1);
    clf;
    plot(log10(rmsei),'k-o');
    ylabel('\fontsize{13} log_{10}(RMSE perf(n))');
    disp(log10(rmsei));
    pause;

end     