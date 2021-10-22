function [TR, TL] = solveSys_3D(RMarkerPts, RMPs, LMarkerPts, LMPs, signR, signL)

RMPs = [RMPs; 0 0 0];
if signR >= 0
    Y = [RMarkerPts [1; 1; 1]]; % ZebrisMoticon
else
    Y = [RMarkerPts([2 1 3],:) [1; 1; 1]]; % ZebrisMoticon 
end
Y = Y';

X = [RMPs; [1, 1, 1]];

TR = Y*pinv(X, eps);

LMPs = [LMPs; 0 0 0];
% Y = [LMarkerPts(:,1:2) [1; 1; 1]]; ZebrisEnvisible
if signL >= 0
    Y = [LMarkerPts([2 1 3],:) [1; 1; 1]]; % ZebrisMoticon        
else
    Y = [LMarkerPts [1; 1; 1]]; % ZebrisMoticon
end
Y = Y';

X = [LMPs; [1, 1, 1]];

% TL = zeros(4);
TL = Y*pinv(X, eps);

end