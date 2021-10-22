% Run_Me: Output -> structures with matrices containing the estimations per
% frame and for each specified subregion
clear; close all
addpath(genpath('./functions'));
addpath('parameter');

%% Choose proband identifier from ["P1", "P2", "P3"]:
Proband = "P1";

%% Choose method from the following (Method) set:
% Method = ["AR", "ARL", "S2S", "S2SL", "SQ", ...
%           "SQL", "BF", "BFL", "NW", "NWL"];
% !!!  Note that, for proband identifier "P3", !!!
% !!!  there is no measurement data            !!!
% !!!  regarding Method "NWL" and "BFL"        !!!
Method = "AR";

%% Insert proband's mass in [kg]
% P1: mass = 73;
% P2: mass = 74;
% P3: mass = 73;
% If the suffix "L" (represents the load) in the Method specifier is chosen, 
% then mass = mass + 14kg. (The mass of the carried load was around 14kg)
mass = 73;

%% eFSM / BE (Bias Estimation)
Estimation = eFSM_BE(char(Method), char(Proband), mass);