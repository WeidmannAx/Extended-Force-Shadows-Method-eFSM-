function Estimation = eFSM_BE(Method, Proband, mass)
%% Set Up
t0 = 2;
Path = strcat('measurements\', Proband, '_', Method, '.mat');

[Force, FSM, ~, Insole, Moticon] = SetUpData(Path, false, mass);

[param, distanceFactor, Bool, ConvHull, Area] = setUp_parameter;
z = 1.0e+03 *[0.002099777557295,   0.000349453107936,   0.000125232027506,   0.005251556352643,   0.000003563034415,   0.016662696532315,   0.000022659109204,   0.000011527893454,   0.002712187924914,   1.250457930355855,   0.312601470155110,   2.105919050930557,   0.074858050263949,   0.860983983053379,   0.087006473064103,   0.476714885637599,   1.183610080428284,   0.854609696335901,   5.912675194195105,   0.004191551391878,   0.000486870951455,   2.145711828376983,   0.933519808381601,   1.792660526606856,   2.421550758140068];

param.filter_C_tilde = diag([1, 1, 1, 1, z(1:9), 1, 1, 1]);
param.filter_meas_variance = diag(z(10:end));
param.filter_inv_meas_variance = inv(param.filter_meas_variance);

%%
FSM.CoP = FSM.Data.centerOfMass;
t = FSM.CoP(3:end,3)./( - Force.F_K(:,3));
FSM.CoP(3:end,:) = FSM.CoP(3:end,:) + t.*(Force.F_K);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tr = FSM.get_to_bottom_transformed_points(Force.F_K);
load('parameter\massDistributionZatsiorsky.mat');
load('parameter\weights.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TR = eye(4); TL = TR;
MSH_left = []; MSH_right = [];
% right orientation
[~,~,yawR]=quat2angle(FSM.Data.orientation(1,18*4+1 : 19*4));
% left orientation
[~,~,yawL]=quat2angle(FSM.Data.orientation(1,22*4+1 : 23*4));

signL = -sign(yawL);
signR = -sign(yawR);

[TR, TL, Parts, FSM.MSH_left, FSM.MSH_right, ~, ~] = getSubregions(FSM, 1, t0, TR, TL, Bool, true, MSH_left, MSH_right, signL, signR, Moticon);

transScale = (0.8469^2)/9.81;
range = 2:1:length(Force.F_K);

GroundContact.CurrentWeightLeft = 1;
GroundContact.CurrentWeightRight = 1;
SumPlate = true;

%% initialize estimation arrays
Estimation.EstimLeft = zeros(length(range), 6); Estimation.EstimRight = Estimation.EstimLeft;
BiasLeft = zeros(length(range), 16); BiasRight = BiasLeft;

FSM.Curve.weightL = weightL;
FSM.Curve.weightR = weightR;
Bool.load_weight_curves = 0;

%% set up for linear regression
% polynomial degree; initial weight parameter and covariances;
polynom_degree = 2;
Sigma_t_L = 1*repmat(eye(polynom_degree+1),16,1);
Sigma_t_R = 1*repmat(eye(polynom_degree+1),16, 1);
w_t_L = 0.1*ones(polynom_degree+1,16);
w_t_R = 0.1*ones(polynom_degree+1,16);

%% main for-loop
for i = range(1):range(50)
    %% 
    [Parts, TR, TL, ConvHull, GroundContact, Estimation, FSM] = Extended_Force_Shadows_Method(i, t0, range, Bool, Moticon, FSM, Force, tr, distribution, param, distanceFactor, Parts, TR, TL, ConvHull, SumPlate, GroundContact, Estimation, transScale);
    %% Bias Estimation (BE)
    if Bool.Bias_Estimation
        %% left foot
        [bias,Sigma_t_L, w_t_L]=Bias_Estimation(Insole, i, t0, Estimation.EstimLeft(i-range(1)+1,:), Area.left, 1, param, Sigma_t_L, w_t_L, polynom_degree);
        % Current estimated bias
        BiasLeft(i-range(1)+1,:) = bias';
        
        %% right foot
        [bias,Sigma_t_R, w_t_R]=Bias_Estimation(Insole, i, t0, Estimation.EstimRight(i-range(1)+1,:), Area.right, 0, param, Sigma_t_R, w_t_R, polynom_degree);
        % Current estimated bias
        BiasRight(i-range(1)+1,:) = bias';
    end
end

%% apply current (latest) regression parameters to all data
for i = range
    %% right side
    Phi = ones(16,polynom_degree+1);
    X = (Insole.pressuresRight(i+t0, :) .* (Area.right / 9.81))';
    for k = 1:polynom_degree
        Phi(:,k+1) = Phi(:,k).*X;
    end
    
    for k = 1:16
        x_BiasRight(k,1) = Phi(k,:)*w_t_R(:,k);
    end
    % additional information about sensor pad activity
    idxnull = find(Insole.pressuresRight(i+t0, :).*Area.right < .4*9.81);
    x_BiasRight(idxnull) = Insole.pressuresRight(i+t0, idxnull)';
    BiasRight(i-range(1)+1,:) = x_BiasRight';
    
    %% left side
    Phi = ones(16,polynom_degree+1);
    X = (Insole.pressuresLeft(i+t0, :) .* (Area.left / 9.81))';
    for k = 1:polynom_degree
        Phi(:,k+1) = Phi(:,k).*X;
    end
    for k = 1:16
        x_BiasLeft(k,1) = Phi(k,:)*w_t_L(:,k);
    end
    % additional information about sensor pad activity
    idxnull = find(Insole.pressuresLeft(i+t0, :).*Area.left < .4*9.81);
    x_BiasLeft(idxnull) = Insole.pressuresLeft(i+t0, idxnull)';
    BiasLeft(i-range(1)+1,:) = x_BiasLeft';
end
%% Insole Correction
if Bool.Bias_Estimation
    Estimation = Shift_By_Bias(range, Insole, t0, Area, BiasLeft, BiasRight, Estimation);
end
Estimation.RawLeft = Estimation.EstimLeft; Estimation.RawRight = Estimation.EstimRight;
Estimation.EstimLeft = filtfilt(designfilt('lowpassiir','FilterOrder',1, 'HalfPowerFrequency',.2,'DesignMethod','butter'), Estimation.EstimLeft);
Estimation.EstimRight = filtfilt(designfilt('lowpassiir','FilterOrder',1, 'HalfPowerFrequency',.2,'DesignMethod','butter'), Estimation.EstimRight);
%% Moticon geometry
Estimation.EstimLeft_Moticon_geometry = Insole.pressuresLeft(range+t0, :).*(Area.left / 9.81) - BiasLeft(1:length(range),:);
Estimation.EstimRight_Moticon_geometry = Insole.pressuresRight(range+t0, :).*(Area.right / 9.81) - BiasRight(1:length(range),:);

end