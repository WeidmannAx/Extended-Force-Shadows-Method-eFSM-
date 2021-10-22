function [f, FSM] = GetCurve(FSM, i, t0, Force, tr, distribution, a, CoP, distanceFactor, Bool, extremities)

[XY, x, y, f] = FSM.getForceShadow(i, t0, Force.F_K, tr, distribution, FSM.CoP, a, CoP, distanceFactor, Bool, extremities);

FSM.Curve.fx = f;
FSM.Curve.x = x;
FSM.Curve.y = y;
FSM.Curve.XY = XY;

end