function [LMarkerPts, RMarkerPts, FtParts] = FootPartitioning(Bool)
%% Foot Mesh
XFtR = [91.14;  96.35;  99.56; 101.96; 103.30; 104.84; 106.51; 106.77; 106.57; 106.57; ...
    106.64; 105.04; 105.24; 103.37;  97.55;  91.34;  87.33; 83.52;  74.77;  72.76;  71.36; ...
    72.83;  69.96;  70.29;  71.76;  74.37;  78.44;  79.04;  81.25;  85.79;  91.14];

YFtR = [73.50;  74.17;  76.71;  82.85;  99.96; 105.77; 111.18; 124.81; 127.89; 138.78; ...
    144.39; 147.47; 156.75; 159.63; 166.51; 170.72; 172.19; 173.26; 169.65; 167.18; 163.90; ...
    156.22; 152.28; 147.07; 142.92; 138.71; 107.38;  81.38;  77.24;  74.50;  73.50];
if Bool.healthy % else flatfoot
    XFtR = [XFtR(1:26);  79.31;  85.66;  89.07;  90.74;  91.34;  89.4;   86.66; 81.98; 80.11; 79.04; XFtR(28:end)];
    YFtR = [YFtR(1:26); 135.51; 132.57; 128.36; 124.35; 116.93; 109.78; 100.49; 92.94; 89.27; 84.06; YFtR(28:end)];
end
FtParts.WholeFootR = [XFtR,YFtR]; FtParts.WholeFootL = [-XFtR,YFtR];
FtParts.Three_Dim_WholeFootR = [FtParts.WholeFootR, zeros(size(FtParts.WholeFootR,1), 1)];
FtParts.Three_Dim_WholeFootL = [FtParts.WholeFootL, zeros(size(FtParts.WholeFootL,1), 1)];

%% Marker points
MPR = [76.17, 103.63, 90.6; 149.6, 141.65, 80.85];

LMarkerPts = [-MPR(1,[2 1 3]);MPR(2,[2 1 3])];
RMarkerPts = MPR;

%% Regions
%% Toe / Antetarsus?
XToeR = [105.04; 105.24; 103.37;  97.55;  91.34;  83.52;  74.77;  72.76;  71.36;  72.83; ...
    86.73;  96.28; 105.04];
YToeR = [147.47; 156.75; 159.63; 166.51; 170.72; 173.26; 169.65; 167.18; 163.90; 156.22; ...
    159.23; 154.68; 147.47];
XToeR_R = [XToeR(1:5);  87.33;  88.4 ; XToeR(12:end)]; YToeR_R = [YToeR(1:5); 172.19; 158.76; YToeR(12:end)];
XToeR_L = [ 87.33; XToeR(6:11);  88.4  ;  87.33]; YToeR_L = [172.19; YToeR(6:11); 158.76 ; 172.19];
% right foot
FtParts.ZebrisToeLR = [XToeR_L,YToeR_L]; FtParts.ZebrisToeRR = [XToeR_R,YToeR_R];
% left foot
FtParts.ZebrisToeLL = [-XToeR_R,YToeR_R]; FtParts.ZebrisToeRL = [-XToeR_L,YToeR_L];

% partition toes into left and right
XToeR = [105.24; 103.37;  97.55;  91.34;  83.52;  74.77;  72.76;  71.36;  72.43;  76.37;  ...
    79.38;  80.65;  82.72;  92.48;  95.62; 101.09; 105.1;  105.24];
YToeR = [156.75; 159.63; 166.51; 170.72; 173.26; 169.65; 167.18; 163.90; 160.43; 159.89;  ...
    160.76; 162.10; 166.24; 164.64; 160.63; 154.62; 154.68; 156.75];
% partition Toes into left and right
XToeR_R = [XToeR(1:4);  87.33;  87.6 ; XToeR(14:end)]; YToeR_R = [YToeR(1:4); 172.19; 165.44; YToeR(14:end)];
XToeR_L = [ 87.33; XToeR(5:13);  87.6  ;  87.33]; YToeR_L = [172.19; YToeR(5:13); 165.44 ; 172.19];
% right foot
FtParts.ToeR = [XToeR,YToeR]; 
FtParts.Three_Dim_ToeR = [FtParts.ToeR, zeros(size(FtParts.ToeR,1), 1)];
FtParts.ToeRR = [XToeR_R,YToeR_R];
FtParts.ToeLR = [XToeR_L,YToeR_L];
% left foot
FtParts.ToeL = [-XToeR,YToeR]; 
FtParts.Three_Dim_ToeL = [FtParts.ToeL, zeros(size(FtParts.ToeL,1), 1)];
FtParts.ToeLL = [-XToeR_R,YToeR_R];
FtParts.ToeRL = [-XToeR_L,YToeR_L];

%% Metatarsus
XMetaR = [106.57; 106.57; 106.64; 105.04;  96.28;  86.73;  72.83;  69.96;  70.29;  71.76;  ...
    74.37;  84.92;  96.28; 106.57];
YMetaR = [127.89; 138.78; 144.39; 147.47; 154.68; 159.23; 156.22; 152.28; 147.07; 142.92; ...
    138.71; 138.04; 133.97; 127.89];
% partition metatarsus into left and right
XMetaR_R = [XMetaR(1:5);  88.40;  89.90; XMetaR(13:end)]; YMetaR_R = [YMetaR(1:5); 158.76; 136.00; YMetaR(13:end)];
XMetaR_L = [ 88.40; XMetaR(6:12);  89.90;  88.40]; YMetaR_L = [158.76; YMetaR(6:12); 136.00; 158.76];

% right foot
FtParts.MetaR = [XMetaR,YMetaR];
FtParts.Three_Dim_MetaR = [FtParts.MetaR, zeros(size(FtParts.MetaR,1), 1)];
FtParts.MetaRR = [XMetaR_R,YMetaR_R];
FtParts.MetaLR = [XMetaR_L,YMetaR_L]; 
% left foot
FtParts.MetaL = [-XMetaR,YMetaR];
FtParts.Three_Dim_MetaL = [FtParts.MetaL, zeros(size(FtParts.MetaL,1), 1)];
FtParts.MetaLL = [-XMetaR_R,YMetaR_R];
FtParts.MetaRL = [-XMetaR_L,YMetaR_L];


%% Mid/Arch
XArchR = [104.84; 106.51; 106.77; 106.57;  96.28;  84.92;  74.37;  78.44;  83.12;  89.40; ...
    96.15; 104.84];
YArchR = [105.77; 111.18; 124.81; 127.89; 133.97; 138.04; 138.71; 107.38; 110.72; 109.78; ...
    106.24; 105.77];
if Bool.healthy
    XArchR = [XArchR(1:7); 79.31; 85.66; 89.07; 90.74; 91.34; XArchR(10:end)];
    YArchR = [YArchR(1:7); 135.51; 132.57; 128.36; 124.35; 116.93; YArchR(10:end)];
end
% right / left foot
FtParts.ArchR = [XArchR, YArchR]; 
FtParts.Three_Dim_ArchR = [FtParts.ArchR, zeros(size(FtParts.ArchR,1), 1)];
FtParts.ArchL = [-XArchR, YArchR];
FtParts.Three_Dim_ArchL = [FtParts.ArchL, zeros(size(FtParts.ArchL,1), 1)];

%% Heel
XHeelR = [91.14;  96.35;  99.56; 101.96; 103.30; 104.84;  96.15;  89.40;  83.12;  78.44; ...
    79.04;  81.25;  85.79;  91.14];

YHeelR = [73.50;  74.17;  76.71;  82.85;  99.96; 105.77; 106.24; 109.78; 110.72; 107.38; ...
    81.38;  77.24;  74.50;  73.50];
if Bool.healthy
    XHeelR = [XHeelR(1:8); 86.66; 81.98; 80.11; 79.04; XHeelR(11:end)];
    YHeelR = [YHeelR(1:8); 100.49; 92.94; 89.27; 84.06; YHeelR(11:end)];
end
% right / left foot
FtParts.HeelR = [XHeelR, YHeelR]; 
FtParts.Three_Dim_HeelR = [FtParts.HeelR, zeros(size(FtParts.HeelR,1), 1)];
FtParts.HeelL = [-XHeelR, YHeelR];
FtParts.Three_Dim_HeelL = [FtParts.HeelL, zeros(size(FtParts.HeelL,1), 1)];

FtParts.LeftMarkerpts.FirstMeta = [-76.17, 149.6];
FtParts.LeftMarkerpts.ScndMeta = [-85.13, 151.74];
FtParts.LeftMarkerpts.ThirdMeta = [-91.34, 150.07];
FtParts.LeftMarkerpts.FourthMeta = [-97.09, 146.93];
FtParts.LeftMarkerpts.FifthMeta = [-103.63, 141.65];
FtParts.LeftMarkerpts.Calcaneus = [-90.6, 80.85];

FtParts.RightMarkerpts.FirstMeta = [76.17, 149.6];
FtParts.RightMarkerpts.ScndMeta = [85.13, 151.74];
FtParts.RightMarkerpts.ThirdMeta = [91.34, 150.07];
FtParts.RightMarkerpts.FourthMeta = [97.09, 146.93];
FtParts.RightMarkerpts.FifthMeta = [103.63, 141.65];
FtParts.RightMarkerpts.Calcaneus = [90.6, 80.85];

FtParts.RightToes = [[76.91, 165.31] ; [83.79, 169.58] ; [91.01, 167.58] ; [96.62, 163.84] ; [102.83, 157.09]];
FtParts.LeftToes = [[-76.91, 165.31] ; [-83.79, 169.58] ; [-91.01, 167.58] ; [-96.62, 163.84] ; [-102.83, 157.09]];
end