function [nL, nR, weightL, weightR, fL, fR, ConvHull, GroundContact, FSM] = gettingArchWeights(Parts, MSH_left, MSH_right, f, CoM, ConvHull, vdt, a, times, GroundContact, i_current, scale, CoP, Bool, FSM)
lHeel = length(MSH_left.Heel.p(:,1));
lArch = length(MSH_left.Arch.p(:,1));
lMetaL = length(MSH_left.MetaL.p(:,1));
lMetaR = length(MSH_left.MetaR.p(:,1));
lToeL = length(MSH_left.ToeL.p(:,1));
lToeR = length(MSH_left.ToeR.p(:,1));

y_total = f([[MSH_left.Heel.p(:,1); MSH_left.Arch.p(:,1); MSH_left.MetaL.p(:,1); MSH_left.MetaR.p(:,1); MSH_left.ToeL.p(:,1); MSH_left.ToeR.p(:,1)]; ...
    [MSH_right.Heel.p(:,1); MSH_right.Arch.p(:,1); MSH_right.MetaL.p(:,1); MSH_right.MetaR.p(:,1); MSH_right.ToeL.p(:,1); MSH_right.ToeR.p(:,1)]], ...
    [[MSH_left.Heel.p(:,2); MSH_left.Arch.p(:,2); MSH_left.MetaL.p(:,2); MSH_left.MetaR.p(:,2); MSH_left.ToeL.p(:,2); MSH_left.ToeR.p(:,2)]; ...
    [MSH_right.Heel.p(:,2); MSH_right.Arch.p(:,2); MSH_right.MetaL.p(:,2); MSH_right.MetaR.p(:,2); MSH_right.ToeL.p(:,2); MSH_right.ToeR.p(:,2)]]);
y = y_total(1 : length(y_total)/2);

fL.Heel = y(1 : lHeel);
fL.Arch = y(lHeel+1 : lHeel+lArch);
fL.MetaL = y(lHeel+lArch+1 : lHeel+lArch+lMetaL);
fL.MetaR = y(lHeel+lArch+lMetaL+1 : lHeel+lArch+lMetaL+lMetaR);
fL.ToeL = y(lHeel+lArch+lMetaL+lMetaR+1 : lHeel+lArch+lMetaL+lMetaR+lToeL);
fL.ToeR = y(lHeel+lArch+lMetaL+lMetaR+lToeL+1 : lHeel+lArch+lMetaL+lMetaR+lToeL+lToeR);

lMetaLR = length(MSH_right.MetaL.p(:,1));
lMetaRR = length(MSH_right.MetaR.p(:,1));
lToeLR = length(MSH_right.ToeL.p(:,1));
lToeRR = length(MSH_right.ToeR.p(:,1));

y = y_total(length(y_total)/2 + 1 : end);

fR.Heel = y(1 : lHeel);
fR.Arch = y(lHeel+1 : lHeel+lArch);
fR.MetaL = y(lHeel+lArch+1 : lHeel+lArch+lMetaLR);
fR.MetaR = y(lHeel+lArch+lMetaLR+1 : lHeel+lArch+lMetaLR+lMetaRR);
fR.ToeL = y(lHeel+lArch+lMetaLR+lMetaRR+1 : lHeel+lArch+lMetaLR+lMetaRR+lToeLR);
fR.ToeR = y(lHeel+lArch+lMetaLR+lMetaRR+lToeLR+1 : lHeel+lArch+lMetaLR+lMetaRR+lToeLR+lToeRR);

%% fitting curve to constructed data
if Bool.Arch_Model
    weightL = FSM.Curve.weightL;
    weightR = FSM.Curve.weightR;
else
    weightL = @(x,y) (ones(size([x,y],1), 1));
    weightR = @(x,y) (ones(size([x,y],1), 1));
end
%% pointwise evaluating weights + integral conservation condition
%% Left Foot
weightLeft = weightL([MSH_left.Heel.p(:,1); MSH_left.Arch.p(:,1); MSH_left.MetaL.p(:,1); MSH_left.MetaR.p(:,1); MSH_left.ToeL.p(:,1);  MSH_left.ToeR.p(:,1)], ...
    [MSH_left.Heel.p(:,2); MSH_left.Arch.p(:,2); MSH_left.MetaL.p(:,2); MSH_left.MetaR.p(:,2); MSH_left.ToeL.p(:,2); MSH_left.ToeR.p(:,2)]);
weightLeft(weightLeft<0) = 0;
valsL = [fL.Heel; fL.Arch; fL.MetaL; fL.MetaR; fL.ToeL; fL.ToeR];
wL = sum(valsL)/sum(valsL .* weightLeft);
%% Right Foot
weightRight = weightR([MSH_right.Heel.p(:,1); MSH_right.Arch.p(:,1); MSH_right.MetaL.p(:,1); MSH_right.MetaR.p(:,1); MSH_right.ToeL.p(:,1); MSH_right.ToeR.p(:,1)], ...
    [MSH_right.Heel.p(:,2); MSH_right.Arch.p(:,2); MSH_right.MetaL.p(:,2); MSH_right.MetaR.p(:,2); MSH_right.ToeL.p(:,2); MSH_right.ToeR.p(:,2)]);
weightRight(weightRight<0) = 0;
valsR = [fR.Heel; fR.Arch; fR.MetaL; fR.MetaR; fR.ToeL; fR.ToeR];
wR = sum(valsR)/sum(valsR .* weightRight);

%% Hip Flexor Model
kL = convhull(Parts.WholeFootL);
kR = convhull(Parts.WholeFootR);

if Bool.CoM_Correction
    CoM = 0.5*(CoM(1:2) + CoP);
end

if Bool.Hip_Flexor_Model
    if a*ConvHull.limit < 1
        VelPredict = CoM + vdt;
    else
        VelPredict = CoM + vdt.*(1 : a*ConvHull.limit)';%2*ConvHull.limit)';
    end
    
    if inpolygon(CoM(1), CoM(2), Parts.WholeFootL(kL,1), Parts.WholeFootL(kL,2)) || sum(inpolygon(VelPredict(:,1), VelPredict(:,2), Parts.WholeFootL(kL,1), Parts.WholeFootL(kL,2))) > 0
        %% if CoM or one of the predicted CoM-positions is contained in the left foot
        if inpolygon(CoM(1), CoM(2), Parts.WholeFootL(kL,1), Parts.WholeFootL(kL,2)) && ~min(inpolygon(VelPredict(:,1), VelPredict(:,2), Parts.WholeFootL(kL,1), Parts.WholeFootL(kL,2)))
            ConvHull.Left = ConvHull.Left - times * (ConvHull.Left > - ConvHull.limit);
            ConvHull.Right = ConvHull.Right - times * (ConvHull.Right > - ConvHull.limit);
        else
            ConvHull.Right = ConvHull.Right + times * (ConvHull.Right < ConvHull.limit);
            ConvHull.Left = ConvHull.Left - times * (ConvHull.Left > - ConvHull.limit);
        end
    elseif inpolygon(CoM(1), CoM(2), Parts.WholeFootR(kR,1), Parts.WholeFootR(kR,2)) || sum(inpolygon(VelPredict(:,1), VelPredict(:,2), Parts.WholeFootR(kR,1), Parts.WholeFootR(kR,2)))>0
        %% if CoM or one of the predicted CoM-positions is contained in the right foot
        if inpolygon(CoM(1), CoM(2), Parts.WholeFootR(kR,1), Parts.WholeFootR(kR,2)) && ~min(inpolygon(VelPredict(:,1), VelPredict(:,2), Parts.WholeFootR(kR,1), Parts.WholeFootR(kR,2)))
            ConvHull.Left = ConvHull.Left - times * (ConvHull.Left > - ConvHull.limit);
            ConvHull.Right = ConvHull.Right - times * (ConvHull.Right > - ConvHull.limit);
        else
            ConvHull.Left = ConvHull.Left + times * (ConvHull.Left < ConvHull.limit);
            ConvHull.Right = ConvHull.Right - times * (ConvHull.Right > - ConvHull.limit);
        end
    else
        ConvHull.Left = ConvHull.Left - times * (ConvHull.Left > - ConvHull.limit);
        ConvHull.Right = ConvHull.Right - times * (ConvHull.Right > - ConvHull.limit);
    end
    
    if abs(ConvHull.Right) > ConvHull.limit
        ConvHull.Right = sign(ConvHull.Right)*ConvHull.limit;
    elseif abs(ConvHull.Left) > ConvHull.limit
        ConvHull.Left = sign(ConvHull.Left)*ConvHull.limit;
    end
    
    sR = -1/(2*ConvHull.limit)*ConvHull.Right + 0.5;
    sL = -1/(2*ConvHull.limit)*ConvHull.Left + 0.5;
    
    nR = wR .* (sR + (1-sL));
    nL = wL .* (sL + (1-sR));
else
    nR = wR;
    nL = wL;
end
ConvHull.sL = sL;
ConvHull.sR = sR;
end