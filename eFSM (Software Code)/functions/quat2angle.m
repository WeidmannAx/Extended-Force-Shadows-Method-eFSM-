function [roll, pitch, yaw] = quat2angle(q)
roll = atan2(2*(q(1)*q(2) + q(3)*q(4)), 1 - 2*(q(2)^2 + q(3)^2));

tmp = 2*(q(1)*q(3) - q(4)*q(2));
if abs(tmp) >= 1
    pitch = pi/2*sign(tmp);
else
    pitch = asin(tmp);
end

yaw = atan2(2*(q(1)*q(4) + q(2)*q(3)), 1 - 2*(q(3)^2 + q(4)^2));
end