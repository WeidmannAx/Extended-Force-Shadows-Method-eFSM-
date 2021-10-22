function [Int] = get_integral_vals(MSH_left, MSH_right, weightL, weightR, ArchBool, nL, nR, fL, fR, Bool_Hip_Flexor_Model)
Int.Heel_L = - compute_integral(MSH_left.Heel,weightL, ArchBool, nL, fL.Heel, Bool_Hip_Flexor_Model);
Int.Heel_R = - compute_integral(MSH_right.Heel,weightR, ArchBool, nR, fR.Heel, Bool_Hip_Flexor_Model);

Int.Arch_L = - compute_integral(MSH_left.Arch,weightL, ArchBool, nL, fL.Arch, Bool_Hip_Flexor_Model);
Int.Arch_R = - compute_integral(MSH_right.Arch,weightR, ArchBool, nR, fR.Arch, Bool_Hip_Flexor_Model);

Int.MetaL_L = - compute_integral(MSH_left.MetaL,weightL, ArchBool, nL, fL.MetaL, Bool_Hip_Flexor_Model);
Int.MetaR_L = - compute_integral(MSH_left.MetaR,weightL, ArchBool, nL, fL.MetaR, Bool_Hip_Flexor_Model);
Int.MetaL_R = - compute_integral(MSH_right.MetaL,weightR, ArchBool, nR, fR.MetaL, Bool_Hip_Flexor_Model);
Int.MetaR_R = - compute_integral(MSH_right.MetaR,weightR, ArchBool, nR, fR.MetaR, Bool_Hip_Flexor_Model);

Int.ToeL_L = - compute_integral(MSH_left.ToeL,weightL, ArchBool, nL, fL.ToeL, Bool_Hip_Flexor_Model);
Int.ToeR_L = - compute_integral(MSH_left.ToeR,weightL, ArchBool, nL, fL.ToeR, Bool_Hip_Flexor_Model);
Int.ToeL_R = - compute_integral(MSH_right.ToeL,weightR, ArchBool, nR, fR.ToeL, Bool_Hip_Flexor_Model);
Int.ToeR_R = - compute_integral(MSH_right.ToeR,weightR, ArchBool, nR, fR.ToeR, Bool_Hip_Flexor_Model);
end