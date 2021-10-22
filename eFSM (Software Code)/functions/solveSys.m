function [TR, TL] = solveSys(RMarkerPts, RMPs, LMarkerPts, LMPs, signR, signL)

% Y = [RMarkerPts([2 1 3],1:2) [1; 1; 1]]; ZebrisEnvisible
if signR >= 0
    Y = [RMarkerPts(:,1:2) [1; 1; 1]]; % ZebrisMoticon
else
    Y = [RMarkerPts([2 1 3],1:2) [1; 1; 1]]; % ZebrisMoticon 
end
Y = Y';

X = [RMPs; [1, 1, 1]];

TR = zeros(4);
TR([1:2, 4],[1:2, 4]) = Y*pinv(X, eps);

% Y = [LMarkerPts(:,1:2) [1; 1; 1]]; ZebrisEnvisible
if signL >= 0
    Y = [LMarkerPts([2 1 3],1:2) [1; 1; 1]]; % ZebrisMoticon        
else
    Y = [LMarkerPts(:,1:2) [1; 1; 1]]; % ZebrisMoticon
end
Y = Y';

X = [LMPs; [1, 1, 1]];

TL = zeros(4);
TL([1:2, 4],[1:2, 4]) = Y*pinv(X, eps);

% disp(strcat('det R = ', num2str(det(TR(1:2,1:2)))));
% disp(strcat('s = ', num2str(TR(4,4))));
% disp(strcat('det R = ', num2str(det(TL(1:2,1:2)))));
% disp(strcat('s = ', num2str(TL(4,4))));

end