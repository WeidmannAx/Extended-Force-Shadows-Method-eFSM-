function Int = IsContactRegion(Int, Parts, Prior_Parts, dt, GroundContact)

MeanHeights = [mean(Parts.Three_Dim_HeelL(:,3)), mean(Parts.Three_Dim_HeelR(:,3)), mean(Parts.Three_Dim_ArchL(:,3)), mean(Parts.Three_Dim_ArchR(:,3)), ...
    mean(Parts.Three_Dim_MetaL(:,3)), mean(Parts.Three_Dim_MetaL(:,3)), mean(Parts.Three_Dim_MetaR(:,3)), mean(Parts.Three_Dim_MetaR(:,3)), ...
    mean(Parts.Three_Dim_ToeL(:,3)), mean(Parts.Three_Dim_ToeL(:,3)), mean(Parts.Three_Dim_ToeR(:,3)), mean(Parts.Three_Dim_ToeR(:,3))];
IsSmall = MeanHeights < 0.03;
%% velocity
PriorHeights = [mean(Prior_Parts.Three_Dim_HeelL(:,3)), mean(Prior_Parts.Three_Dim_HeelR(:,3)), mean(Prior_Parts.Three_Dim_ArchL(:,3)), mean(Prior_Parts.Three_Dim_ArchR(:,3)), ...
    mean(Prior_Parts.Three_Dim_MetaL(:,3)), mean(Prior_Parts.Three_Dim_MetaL(:,3)), mean(Prior_Parts.Three_Dim_MetaR(:,3)), mean(Prior_Parts.Three_Dim_MetaR(:,3)), ...
    mean(Prior_Parts.Three_Dim_ToeL(:,3)), mean(Prior_Parts.Three_Dim_ToeL(:,3)), mean(Prior_Parts.Three_Dim_ToeR(:,3)), mean(Prior_Parts.Three_Dim_ToeR(:,3))];
VelocityRegions = (MeanHeights - PriorHeights)/dt;
IsSlow = VelocityRegions < 0.8;

IsBoth = IsSmall & IsSlow;

Int.Heel_L = Int.Heel_L * IsBoth(1);
Int.Heel_R = Int.Heel_R * IsBoth(2);
Int.Arch_L = Int.Arch_L * IsBoth(3);
Int.Arch_R = Int.Arch_R * IsBoth(4);
Int.MetaL_L = Int.MetaL_L * IsBoth(5);
Int.MetaR_L = Int.MetaR_L * IsBoth(6);
Int.MetaL_R = Int.MetaL_R * IsBoth(7);
Int.MetaR_R = Int.MetaR_R * IsBoth(8);
Int.ToeL_L = Int.ToeL_L * IsBoth(9);
Int.ToeR_L = Int.ToeR_L * IsBoth(10);
Int.ToeL_R = Int.ToeL_R * IsBoth(11);
Int.ToeR_R = Int.ToeR_R * IsBoth(12);


%% Before, the next part was located at the end of "get_integral_vals" 
Weight = Int.Heel_L + Int.Heel_R + ...
    Int.Arch_L + Int.Arch_R + ...
    Int.MetaL_L + Int.MetaL_R + Int.MetaR_L + Int.MetaR_R + ...
    Int.ToeL_L + Int.ToeL_R + Int.ToeR_L + Int.ToeR_R;

% TotalL = Int.Heel_L + Int.Arch_L + Int.MetaL_L + Int.MetaR_L + Int.ToeL_L + Int.ToeR_L; TotalR = Int.Heel_R + Int.Arch_R + Int.MetaL_R + Int.MetaR_R + Int.ToeL_R + Int.ToeR_R;
Int.Heel_L = Int.Heel_L/Weight*GroundContact.CurrentWeightLeft; 
Int.Heel_R = Int.Heel_R/Weight*GroundContact.CurrentWeightRight;  
Int.Arch_L = Int.Arch_L/Weight*GroundContact.CurrentWeightLeft; 
Int.Arch_R = Int.Arch_R/Weight*GroundContact.CurrentWeightRight;
Int.MetaL_L = Int.MetaL_L/Weight*GroundContact.CurrentWeightLeft; 
Int.MetaR_L = Int.MetaR_L/Weight*GroundContact.CurrentWeightLeft; 

Int.MetaL_R = Int.MetaL_R/Weight*GroundContact.CurrentWeightRight; 
Int.MetaR_R = Int.MetaR_R/Weight*GroundContact.CurrentWeightRight;
Int.ToeL_L = Int.ToeL_L/Weight*GroundContact.CurrentWeightLeft; 
Int.ToeR_L = Int.ToeR_L/Weight*GroundContact.CurrentWeightLeft; 
Int.ToeL_R = Int.ToeL_R/Weight*GroundContact.CurrentWeightRight; 
Int.ToeR_R = Int.ToeR_R/Weight*GroundContact.CurrentWeightRight;

