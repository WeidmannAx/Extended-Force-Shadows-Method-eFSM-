function [funcL, funcR, fitresultR] = createSurfaceFit(Parts, vis, scale)
%CREATEFIT(XFIT,YFIT,ZFIT)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : Xfit
%      Y Input : Yfit
%      Z Output: Zfit
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 02-Jul-2019 07:17:00

% scale = 0.1633;
% heightL = norm(Parts.LeftMarkerpts.Calcaneus-Parts.LeftMarkerpts.FifthMeta)/scale;
% heightR = norm(Parts.RightMarkerpts.Calcaneus-Parts.RightMarkerpts.FifthMeta)/scale;
heightL = scale; 
heightR = scale;

FitVariablesL =[ [Parts.LeftMarkerpts.FirstMeta, 1]; [Parts.LeftMarkerpts.ScndMeta, 1]; [Parts.LeftMarkerpts.ThirdMeta, 1]; ...
    [Parts.LeftMarkerpts.FourthMeta, 1]; [Parts.LeftMarkerpts.FifthMeta, 1]; [Parts.LeftMarkerpts.Calcaneus, 1]; ...
    [(Parts.LeftMarkerpts.Calcaneus + Parts.LeftMarkerpts.FifthMeta)/2, heightL]; ...
    [(Parts.LeftMarkerpts.Calcaneus + Parts.LeftMarkerpts.FirstMeta)/2, heightL]]; %...
%     [Parts.LeftToes(1:5,:), linspace(1,.9,5)']];
%                 [Parts.LeftToes, ones(length(Parts.LeftToes),1)]];
FitVariablesR =[ [Parts.RightMarkerpts.FirstMeta, 1]; [Parts.RightMarkerpts.ScndMeta, 1]; [Parts.RightMarkerpts.ThirdMeta, 1]; ...
    [Parts.RightMarkerpts.FourthMeta, 1]; [Parts.RightMarkerpts.FifthMeta, 1]; [Parts.RightMarkerpts.Calcaneus, 1]; ...
    [(Parts.RightMarkerpts.Calcaneus + Parts.RightMarkerpts.FifthMeta)/2, heightR]; ...
    [(Parts.RightMarkerpts.Calcaneus + Parts.RightMarkerpts.FirstMeta)/2, heightR]]; ...
%     [Parts.RightToes(1:5,:), linspace(1,.9,5)']]; % vorher 0.25

XfitL = FitVariablesL(:,1);
YfitL = FitVariablesL(:,2);
ZfitL = FitVariablesL(:,3);

XfitR = FitVariablesR(:,1);
YfitR = FitVariablesR(:,2);
ZfitR = FitVariablesR(:,3);

%% Fit: 'untitled fit 1'.
[xDataL, yDataL, zDataL] = prepareSurfaceData( XfitL, YfitL, ZfitL );

% Set up fittype and options.
ft = fittype( 'poly22' );
% ft = fittype( 'poly23' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Normalize = 'on';
% opts.Robust = 'LAR';

% Fit model to data.
% [fitresultL, gof] = fit( [xDataL, yDataL], zDataL, ft );
[fitresultL, gof] = fit( [xDataL, yDataL], zDataL, ft, opts );

% figure
funcL = @(x,y) fitresultL(x,y);

%% Fit: 'untitled fit 1'.
[xDataR, yDataR, zDataR] = prepareSurfaceData( XfitR, YfitR, ZfitR );

[fitresultR, gof] = fit( [xDataR, yDataR], zDataR, ft, opts );
funcR = @(x,y) fitresultR(x,y);

if vis
    % Plot fit with data.
    figure(2);
    subplot(121)
    hold on; plot(Parts.WholeFootL(:,1),Parts.WholeFootL(:,2),'k')
    plot( fitresultL, [xDataL, yDataL], zDataL );
    alpha 0.25
    shading interp
    subplot(122)
    hold on; plot(Parts.WholeFootR(:,1),Parts.WholeFootR(:,2),'k')
    plot( fitresultR, [xDataR, yDataR], zDataR );
    legend('Zfit vs. Xfit, Yfit', 'Location', 'NorthEast' );
    % Label axes
    xlabel Xfit
    ylabel Yfit
    zlabel Zfit
    grid on
    view( -98.3, 26.0 );
    alpha 0.25
    shading interp
    figure(1);
end

