function [TR, TL, Parts, MSH_left, MSH_right, LMarkerPts, RMarkerPts] = getSubregions(MVN, i, t0, TR, TL, Bool, getTriang, MSH_left_prior, MSH_right_prior, signL, signR, Moticon)
FootR = [52:54; 55:57];
FootL = [64:66; 67:69];

try
    [LMarkerPts, RMarkerPts] = MVN.getMarkerPoints(i+t0);
catch
    [LMarkerPts] = getFeet(MVN, FootL, i+t0, 1, 0); % bis 11.12.2019 hatte ich i+t0 = i
    [RMarkerPts] = getFeet(MVN, FootR, i+t0, -1, 1);% bis 11.12.2019 hatte ich i+t0 = i
end

if Bool.UseMoticonGeometry
    [LMPs, RMPs, Parts] = FootPartitioning_Moticon(Moticon);
    [LMarkerPts, RMarkerPts, ~] = FootPartitioning(Bool);
    LMarkerPts(:,[1,2]) = LMarkerPts(:,[2,1]);
    LMarkerPts = [LMarkerPts; 0, 0, 0];
    RMarkerPts = [RMarkerPts; 0, 0, 0];
    [Parts.TR, Parts.TL] = solveSys_3D(RMarkerPts', RMPs, LMarkerPts', LMPs, signR, signL);
    MSH_left = []; MSH_right = [];
    TR = []; TL = [];
else
    [LMPs, RMPs, Parts] = FootPartitioning(Bool);

    if getTriang
            load('MeshData.mat', 'MSH_left', 'MSH_right')
    else
        
        [TR, TL] = solveSys_3D(RMarkerPts, RMPs, LMarkerPts, LMPs, signR, signL);
        
        
        labels = ["Heel", "Arch", "MetaL", "MetaR", "ToeL", "ToeR", "WholeFoot"];
        %     for j = 1:length(labels)
        [MSH_left.Heel, MSH_right.Heel] = transform_points(TL, TR, MSH_left_prior.Heel, MSH_right_prior.Heel);
        [MSH_left.Arch, MSH_right.Arch] = transform_points(TL, TR, MSH_left_prior.Arch, MSH_right_prior.Arch);
        [MSH_left.MetaL, MSH_right.MetaL] = transform_points(TL, TR, MSH_left_prior.MetaL, MSH_right_prior.MetaL);
        [MSH_left.MetaR, MSH_right.MetaR] = transform_points(TL, TR, MSH_left_prior.MetaR, MSH_right_prior.MetaR);
        [MSH_left.ToeL, MSH_right.ToeL] = transform_points(TL, TR, MSH_left_prior.ToeL, MSH_right_prior.ToeL);
        [MSH_left.ToeR, MSH_right.ToeR] = transform_points(TL, TR, MSH_left_prior.ToeR, MSH_right_prior.ToeR);
        Parts.Three_Dim_MetaL = (TL(4,4)*TL(1:3,1:3)*Parts.Three_Dim_MetaL' + TL(1:3,4))';
        Parts.Three_Dim_MetaR = (TR(4,4)*TR(1:3,1:3)*Parts.Three_Dim_MetaR' + TR(1:3,4))';
        Parts.Three_Dim_ToeL = (TL(4,4)*TL(1:3,1:3)*Parts.Three_Dim_ToeL' + TL(1:3,4))';
        Parts.Three_Dim_ToeR = (TR(4,4)*TR(1:3,1:3)*Parts.Three_Dim_ToeR' + TR(1:3,4))';
        Parts.Three_Dim_HeelL = (TL(4,4)*TL(1:3,1:3)*Parts.Three_Dim_HeelL' + TL(1:3,4))';
        Parts.Three_Dim_ArchL = (TL(4,4)*TL(1:3,1:3)*Parts.Three_Dim_ArchL' + TL(1:3,4))';
        Parts.Three_Dim_WholeFootL = (TL(4,4)*TL(1:3,1:3)*Parts.Three_Dim_WholeFootL' + TL(1:3,4))';
        Parts.Three_Dim_HeelR = (TR(4,4)*TR(1:3,1:3)*Parts.Three_Dim_HeelR' + TR(1:3,4))';
        Parts.Three_Dim_ArchR = (TR(4,4)*TR(1:3,1:3)*Parts.Three_Dim_ArchR' + TR(1:3,4))';
        Parts.Three_Dim_WholeFootR = (TR(4,4)*TR(1:3,1:3)*Parts.Three_Dim_WholeFootR' + TR(1:3,4))';
        %         end
        Parts.HeelL = (TL(4,4)*TL(1:2,1:2)*Parts.HeelL' + TL(1:2,4))';
        Parts.ArchL = (TL(4,4)*TL(1:2,1:2)*Parts.ArchL' + TL(1:2,4))';
        Parts.MetaLL = (TL(4,4)*TL(1:2,1:2)*Parts.MetaLL' + TL(1:2,4))';
        Parts.MetaRL = (TL(4,4)*TL(1:2,1:2)*Parts.MetaRL' + TL(1:2,4))';
        Parts.ToeLL = (TL(4,4)*TL(1:2,1:2)*Parts.ToeLL' + TL(1:2,4))';
        Parts.ToeRL = (TL(4,4)*TL(1:2,1:2)*Parts.ToeRL' + TL(1:2,4))';
        Parts.WholeFootL = (TL(4,4)*TL(1:2,1:2)*Parts.WholeFootL' + TL(1:2,4))';
        
        Parts.HeelR = (TR(4,4)*TR(1:2,1:2)*Parts.HeelR' + TR(1:2,4))';
        Parts.ArchR = (TR(4,4)*TR(1:2,1:2)*Parts.ArchR' + TR(1:2,4))';
        Parts.MetaLR = (TR(4,4)*TR(1:2,1:2)*Parts.MetaLR' + TR(1:2,4))';
        Parts.MetaRR = (TR(4,4)*TR(1:2,1:2)*Parts.MetaRR' + TR(1:2,4))';
        Parts.ToeLR = (TR(4,4)*TR(1:2,1:2)*Parts.ToeLR' + TR(1:2,4))';
        Parts.ToeRR = (TR(4,4)*TR(1:2,1:2)*Parts.ToeRR' + TR(1:2,4))';
        Parts.WholeFootR = (TR(4,4)*TR(1:2,1:2)*Parts.WholeFootR' + TR(1:2,4))';
        
        Parts.RightMarkerpts.FirstMeta = (TR(4,4)*TR(1:2,1:2)*Parts.RightMarkerpts.FirstMeta' + TR(1:2,4))';
        Parts.RightMarkerpts.ScndMeta = (TR(4,4)*TR(1:2,1:2)*Parts.RightMarkerpts.ScndMeta' + TR(1:2,4))';
        Parts.RightMarkerpts.ThirdMeta = (TR(4,4)*TR(1:2,1:2)*Parts.RightMarkerpts.ThirdMeta' + TR(1:2,4))';
        Parts.RightMarkerpts.FourthMeta = (TR(4,4)*TR(1:2,1:2)*Parts.RightMarkerpts.FourthMeta' + TR(1:2,4))';
        Parts.RightMarkerpts.FifthMeta = (TR(4,4)*TR(1:2,1:2)*Parts.RightMarkerpts.FifthMeta' + TR(1:2,4))';
        Parts.RightMarkerpts.Calcaneus = (TR(4,4)*TR(1:2,1:2)*Parts.RightMarkerpts.Calcaneus' + TR(1:2,4))';
        
        Parts.LeftMarkerpts.FirstMeta = (TL(4,4)*TL(1:2,1:2)*Parts.LeftMarkerpts.FirstMeta' + TL(1:2,4))';
        Parts.LeftMarkerpts.ScndMeta = (TL(4,4)*TL(1:2,1:2)*Parts.LeftMarkerpts.ScndMeta' + TL(1:2,4))';
        Parts.LeftMarkerpts.ThirdMeta = (TL(4,4)*TL(1:2,1:2)*Parts.LeftMarkerpts.ThirdMeta' + TL(1:2,4))';
        Parts.LeftMarkerpts.FourthMeta = (TL(4,4)*TL(1:2,1:2)*Parts.LeftMarkerpts.FourthMeta' + TL(1:2,4))';
        Parts.LeftMarkerpts.FifthMeta = (TL(4,4)*TL(1:2,1:2)*Parts.LeftMarkerpts.FifthMeta' + TL(1:2,4))';
        Parts.LeftMarkerpts.Calcaneus = (TL(4,4)*TL(1:2,1:2)*Parts.LeftMarkerpts.Calcaneus' + TL(1:2,4))';
        
        Parts.RightToes = (TR(4,4)*TR(1:2,1:2)*Parts.RightToes' + TR(1:2,4))';
        
        Parts.LeftToes = (TL(4,4)*TL(1:2,1:2)*Parts.LeftToes' + TL(1:2,4))';
        
    end
end
end

function [points_left, points_right] = transform_points(TL, TR, points_left_prior, points_right_prior)
points_left.p = (TL(4,4)*TL(1:2,1:2)*points_left_prior.p' + TL(1:2,4))';
points_right.p = (TR(4,4)*TR(1:2,1:2)*points_right_prior.p' + TR(1:2,4))';
points_left.bp = (TL(4,4)*TL(1:2,1:2)*points_left_prior.bp' + TL(1:2,4))';
points_right.bp = (TR(4,4)*TR(1:2,1:2)*points_right_prior.bp' + TR(1:2,4))';
points_left.TR = points_left_prior.TR;
points_right.TR = points_right_prior.TR;
end