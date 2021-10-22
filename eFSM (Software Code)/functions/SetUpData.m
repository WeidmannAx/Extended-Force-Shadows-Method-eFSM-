function [Force, MVN, Zeb, Insole, Moticon] = SetUpData(Path, vis, mass)
%% load measurement file
load(Path);

%% Define Moticon object
Insole = Data.Moticon; Moticon = MOTICON; Moticon.data = Insole;
%% Defining AwindaMVN object
MVN = AwindaMVN; MVN.Frames = Data.frames;
MVN.BodyWeight = mass; MVN.Frequency = Data.sampleRate; MVN.Data = Data.MVNX;

F_K = MVN.getForce; 
FNorm = vecnorm(F_K')';
Kinematic_Only = median(nonzeros(FNorm))/9.81;
        
%% Defining Zebris object
Zeb = ZebrisFDM15; Zeb.Frames = Data.frames; Zeb.Frequency = Data.sampleRate; Zeb.Data = Data.Zebris.pressureMatrices;
Zebris = Data.Zebris;

F_Load = Zeb.getForce;
F_Load = F_Load(1:length(MVN.Data.position));
F_Load = F_Load(1:length(F_K) + 2);
if vis
    figure
    plot(F_Load); hold on; title('Pressure Plate Data (Zebris FDM 1.5)');
end
[idxI, ~, ~] = find(F_Load);
F_Load(idxI) = F_Load(idxI) - Kinematic_Only*9.81; F_Load = F_Load.*(F_Load >= 0);
if vis
    plot(F_Load); hold off; legend(strcat('load only (mean: ', num2str(mean(nonzeros(F_Load./9.81))), 'kg)'));
end
if ~isempty(nonzeros(F_Load./9.81))
    Load_Only = median(nonzeros(F_Load./9.81));
else
    Load_Only = 0;
end
[phi, psi, gamma] = MVN.getForceAngles; Load_Vector = [cos(phi), cos(psi), cos(gamma)].*repmat(F_Load(3:end), 1, 3);
Load_Vector = Load_Vector./9.81;
if vis
    figure; subplot(131); plot(Load_Vector(:,1)); subplot(132); plot(Load_Vector(:,2)); subplot(133); plot(Load_Vector(:,3));
end

%% draw Force resulting from kinematics

Force.F_T = mass; Force.F_K = F_K; Force.F_Load = F_Load;
Force.Load_Vector = Load_Vector; Force.FNorm = FNorm;
Force.Load = Load_Only;
end