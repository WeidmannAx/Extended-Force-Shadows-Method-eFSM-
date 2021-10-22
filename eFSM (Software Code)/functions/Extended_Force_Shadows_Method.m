function [Parts, TR, TL, ConvHull, GroundContact, Estimation, FSM] = Extended_Force_Shadows_Method(i, t0, range, Bool, Moticon, FSM, Force, tr, distribution, param, distanceFactor, Parts, TR, TL, ConvHull, SumPlate, GroundContact, Estimation, transScale)
% right orientation
[~,~,yawR]=quat2angle(FSM.Data.orientation(i+t0,18*4+1 : 19*4));
% left orientation
[~,~,yawL]=quat2angle(FSM.Data.orientation(i+t0,22*4+1 : 23*4));
signL = -sign(yawL);
signR = -sign(yawR);

%% CoM correction
if Bool.CoM_Correction
    CoP = totalInsoleCoP(Moticon, FSM, i, t0);
    CoP_1 = totalInsoleCoP(Moticon, FSM, i-1, t0);
else
    CoP = [];
    CoP_1 = [];
end
%% Kinematic / 1st Subplot
[f, FSM] = GetCurve(FSM, i, t0, Force, tr, distribution, param.segment_covariances, CoP, distanceFactor, Bool, param.extremities);

Prior_Parts = Parts;
[TR, TL, Parts, MSH_left, MSH_right, LMarkerPts, RMarkerPts] = getSubregions(FSM, i, t0, TR, TL, Bool, false, FSM.MSH_left, FSM.MSH_right, signL, signR, Moticon);

%% Weighting and Integration
if Bool.CoM_Correction
    v = FSM.Frequency * ((FSM.CoP(i+t0,1:2) + distanceFactor*(FSM.CoP(i+t0,1:2) - CoP(1:2))) - (FSM.CoP(i+t0-1,1:2) + distanceFactor*(FSM.CoP(i+t0-1,1:2) - CoP_1(1:2))));
    COM = FSM.CoP(i+t0,1:2) + distanceFactor*(FSM.CoP(i+t0,1:2) - CoP(1:2));
else
    v = FSM.Frequency * (FSM.CoP(i+t0,1:2) - FSM.CoP(i+t0-1,1:2));
    COM = FSM.CoP(i+t0,1:2);
end

[nL, nR, weightL, weightR, fL, fR, ConvHull, GroundContact, FSM] = gettingArchWeights(Parts, MSH_left, MSH_right, f, COM, ConvHull, 1/FSM.Frequency*v, param.ramp(2), param.ramp(3), GroundContact, i, param.arch_height, CoP, Bool, FSM);

[Int] = get_integral_vals(MSH_left, MSH_right, weightL, weightR, Bool.Arch_Model, nL, nR, fL, fR, Bool.Hip_Flexor_Model);
Int = IsContactRegion(Int, Parts, Prior_Parts, 1/FSM.Frequency, GroundContact);
%% write vector to data
Estimation.EstimLeft(i-range(1)+1,:) = [Int.Heel_L, Int.Arch_L, Int.MetaL_L, Int.MetaR_L, Int.ToeL_L, Int.ToeR_L]*(abs(Force.F_K(i,3)/9.81) + Force.Load);
Estimation.EstimRight(i-range(1)+1,:) = [Int.Heel_R, Int.Arch_R, Int.MetaL_R, Int.MetaR_R, Int.ToeL_R, Int.ToeR_R]*(abs(Force.F_K(i,3)/9.81) + Force.Load);
end

function C = totalInsoleCoP(Moticon, FSM, i, t0)
% right orientation
[~,~,yawR]=quat2angle(FSM.Data.orientation(i+t0,18*4+1 : 19*4));
% left orientation
[~,~,yawL]=quat2angle(FSM.Data.orientation(i+t0,22*4+1 : 23*4));
signL = -sign(yawL);
signR = -sign(yawR);

[LM, RM] = Moticon.getLandmarks_CoP_Based;
LM = LM(:,[2 1])'; RM = [-RM(:,2), RM(:,1)]';
[LMarkerPts, RMarkerPts] = FSM.getMarkerPoints(i+t0);
[TransformCoPR, TransformCoPL] = solveSys(RMarkerPts, RM, LMarkerPts, LM, signR, signL);
CoPL = Moticon.data.CoPLeft(i+t0,[2 1]);
CoPR = [ - Moticon.data.CoPRight(i+t0, 2), Moticon.data.CoPRight(i+t0, 1)];

if CoPR(1) ~= 0 || CoPR(2) ~= 0
    CoPR = (TransformCoPR(4,4)*TransformCoPR(1:2,1:2) * CoPR' + TransformCoPR(1:2,4))';
end
if CoPL(1) ~= 0 || CoPL(2) ~= 0
    CoPL = (TransformCoPL(4,4)*TransformCoPL(1:2,1:2) * CoPL' + TransformCoPL(1:2,4))';
end
PL = sum(Moticon.data.pressuresLeft(i+t0,:));
PR = sum(Moticon.data.pressuresRight(i+t0,:));
C = (PL*CoPL + PR*CoPR)/(PL + PR);
end