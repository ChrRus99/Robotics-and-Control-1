function Rz = R_z(theta)
% Function that returns the elementary rotation matrix around z-axis
Rz = [[cos(theta), -sin(theta), 0];...
      [sin(theta),  cos(theta), 0];...
      [0,           0,          1]];