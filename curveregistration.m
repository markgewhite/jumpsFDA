% ************************************************************************
% Function: curveregistration
% Purpose:  Perform landmark or continuous registration
%
% Parameters:
%       t: timespan (array of time points)
%       xFd: curve to be registered (functional data object)
%       lm: (optional) structure of .mean and .cases  of landmark times
%       basis: number of bases (LM only)
%       order: basis order (LM only)
%       lambda: warping roughness penalty (LM only)
%       xlambda: curve roughness penalty (LM only)
%
% Outputs:
%       xFdReg: registered curves
%       warpFdReg: time warp curves
%       wFdReg: W function warp curves
%
% ************************************************************************


function [ XFdReg, warpFdReg, WFdReg ] = curveregistration( ...
                                            t, XFd, lm, ...
                                            basis, order, ...
                                            lambda, xlambda)

monotonic = true; % time does not go backwards

%  Set up a simply monomial basis for landmark registration
%  This will compute warping functions that interpolate the landmark times

if isempty(lm)
    wBasisReg = create_bspline_basis( [t(1),t(end)], basis, order );
else
    wBasisReg = create_bspline_basis( [t(1),t(end)], basis, order, ...
                                        [ t(1) lm.mean t(end) ]);
end
warpFdReg = fd( zeros(basis,1), wBasisReg );
wFdParReg = fdPar( warpFdReg, 1, lambda );

if isempty(lm)
    % perform continuous registration
    XMeanFd = mean(XFd);
    [ XFdReg, warpFdReg, WFdReg ] = register_fd( XMeanFd, XFd, wFdParReg) ;
    
else
    % perform landmark registration
    [XFdReg,warpFdReg,WFdReg] = landmarkreg( XFd, lm.case, lm.mean, ...
                                wFdParReg,monotonic,xlambda);
                            
end

if ~isempty(plotcases)
    % plot a selection of cases
    % generate the derivatives
    dxFd = deriv(XFd,1);
    dxFdReg = deriv(XFdReg,1);

    % evaluate the functions
    xPts = eval_fd(t,XFd);
    xPtsReg = eval_fd(t,XFdReg);
    dxPts = eval_fd(t,dxFd);
    dxPtsReg = eval_fd(t,dxFdReg);
    warpPtsReg = eval_fd(t,warpFdReg);

    %  plot registered accelerations along with warping functions
    for i = plotcases
        
        clf;
        % plot the curve x and its registered curve
        subplot(3,1,1);
        hold on;
        plot(t,xPtsReg(:,i),'b-', ...
             t,xPts(:,i),'b--',...
             [1,t(end)],[0,0],'b:');
        ymax = max(xPtsReg(:,i),xPts(:,i));
        ymin = min(xPtsReg(:,i),xPts(:,i));
        for j = 1:length(lm.mean)
            plot([lm.mean(j),lm.mean(j)],[ymin,ymax],'b:');
        end
        hold off;
        xlabel('Time t');
        ylabel('Registered x');
        title(['Jump ',num2str(i)]);
        
        % plot the curve dx and its registered curve
        subplot(3,1,2);
        hold on;
        plot(t,dxPtsReg(:,i),'b-', ...
             t,dxPts(:,i),'b--',...
             [1,t(end)],[0,0],'b:');
        ymax = max(dxPtsReg(:,i),dxPts(:,i));
        ymin = min(dxPtsReg(:,i),dxPts(:,i));
        for j = 1:length(lm.mean)
            plot([lm.mean(j),lm.mean(j)],[ymin,ymax],'b:');
        end
        hold off;
        xlabel('Time t');
        ylabel('Registered dx');
        title(['Jump ',num2str(i)]);
        
        % plot the time warp curve
        subplot(3,1,3);
        plot(t,warpPtsReg(:,i),'b-', ...
             [t(1),t(end)],[t(1),t(end)],'b--');
        hold on;
        for j = 1:length(lm.mean)
            plot(lm.mean(j),lm.case(i,j),'o');
        end
        hold off;
        xlabel('Time t');
        ylabel('Registered Time');
        title(['Jump ',num2str(i)]);
        
        pause;
    end
end

end
