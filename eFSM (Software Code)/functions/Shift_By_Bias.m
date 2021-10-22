function [Estimation] = Shift_By_Bias(range, Insole, t0, Area, BiasLeft, BiasRight, Estimation)
p = 0.4400;
p1 = 0.3469;
%% left foot
Insole_Heel = sum((Insole.pressuresLeft(range+t0, 1:6).*(Area.left(1:6) / 9.81) - BiasLeft(1:length(range),1:6)).*[1, 1, 1, 1, 0.3465, 0.1636], 2);
Insole_Arch = sum((Insole.pressuresLeft(range+t0, 5:8).*(Area.left(5:8) / 9.81) - BiasLeft(1:length(range),5:8)).*[0.6535, 0.8364, 0.7716, 0.6438], 2);
Insole_Meta = sum((Insole.pressuresLeft(range+t0, 7:13).*(Area.left(7:13) / 9.81) - BiasLeft(1:length(range),7:13)).*[0.2284, 0.3562, 0.6508, 0.7621, 0.726, 0.5974, 0.4823], 2);
Insole_Toes = sum((Insole.pressuresLeft(range+t0, 9:16).*(Area.left(9:16) / 9.81) - BiasLeft(1:length(range),9:16)).*[0.3492, 0.2379, 0.274, 0.4026, 0.5177, 1, 1, 1], 2);
Estimation.EstimLeft(1:length(range),:) = [Insole_Heel, Insole_Arch, (1-p)*Insole_Meta, p*Insole_Meta, (1-p1)*Insole_Toes, p1*Insole_Toes];

% Estimation.EstimLeft(1:length(range),1) = sum((Insole.pressuresLeft(range+t0,1:4).*(Area.left(1:4)/9.81) - BiasLeft(1:length(range),1:4)), 2);
% Estimation.EstimLeft(1:length(range),2) = sum((Insole.pressuresLeft(range+t0,5:8).*(Area.left(5:8)/9.81) - BiasLeft(1:length(range),5:8)), 2);
% Estimation.EstimLeft(1:length(range),4) = sum((Insole.pressuresLeft(range+t0,9:10).*(Area.left(9:10)/9.81) - BiasLeft(1:length(range),9:10)), 2);
% Estimation.EstimLeft(1:length(range),3) = sum((Insole.pressuresLeft(range+t0,11:13).*(Area.left(11:13)/9.81) - BiasLeft(1:length(range),11:13)), 2);
% Estimation.EstimLeft(1:length(range),6) = sum((Insole.pressuresLeft(range+t0,14).*(Area.left(14)/9.81) - BiasLeft(1:length(range),14)), 2);
% Estimation.EstimLeft(1:length(range),5) = sum((Insole.pressuresLeft(range+t0,15:16).*(Area.left(15:16)/9.81) - BiasLeft(1:length(range),15:16)), 2);

%% right foot
Insole_Heel = sum((Insole.pressuresRight(range+t0, 1:6).*(Area.right(1:6) / 9.81) - BiasRight(1:length(range),1:6)).*[1, 1, 1, 1, 0.3465, 0.1636], 2);
Insole_Arch = sum((Insole.pressuresRight(range+t0, 5:8).*(Area.right(5:8) / 9.81) - BiasRight(1:length(range),5:8)).*[0.6535, 0.8364, 0.7716, 0.6438], 2);
Insole_Meta = sum((Insole.pressuresRight(range+t0, 7:13).*(Area.right(7:13) / 9.81) - BiasRight(1:length(range),7:13)).*[0.2284, 0.3562, 0.6508, 0.7621, 0.726, 0.5974, 0.4823], 2);
Insole_Toes = sum((Insole.pressuresRight(range+t0, 9:16).*(Area.right(9:16) / 9.81) - BiasRight(1:length(range),9:16)).*[0.3492, 0.2379, 0.274, 0.4026, 0.5177, 1, 1, 1], 2);
Estimation.EstimRight(1:length(range),:) = [Insole_Heel, Insole_Arch, p*Insole_Meta, (1-p)*Insole_Meta, p1*Insole_Toes, (1-p1)*Insole_Toes];
% Estimation.EstimRight(1:length(range),1) = sum((Insole.pressuresRight(range+t0,1:4).*(Area.right(1:4)/9.81) - BiasRight(1:length(range),1:4)), 2);
% Estimation.EstimRight(1:length(range),2) = sum((Insole.pressuresRight(range+t0,5:8).*(Area.right(5:8)/9.81) - BiasRight(1:length(range),5:8)), 2);
% Estimation.EstimRight(1:length(range),3) = sum((Insole.pressuresRight(range+t0,9:10).*(Area.right(9:10)/9.81) - BiasRight(1:length(range),9:10)), 2);
% Estimation.EstimRight(1:length(range),4) = sum((Insole.pressuresRight(range+t0,11:13).*(Area.right(11:13)/9.81) - BiasRight(1:length(range),11:13)), 2);
% Estimation.EstimRight(1:length(range),5) = sum((Insole.pressuresRight(range+t0,14).*(Area.right(14)/9.81) - BiasRight(1:length(range),14)), 2);
% Estimation.EstimRight(1:length(range),6) = sum((Insole.pressuresRight(range+t0,15:16).*(Area.right(15:16)/9.81) - BiasRight(1:length(range),15:16)), 2);
end