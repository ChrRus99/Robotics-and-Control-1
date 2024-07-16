% Define the forward kinematics function
function [x, y, phi] = forward_kinematics(theta1, theta2, theta3)
    % This function returns the end-effector position (x, y, phi) based on
    % the joint angles (theta1, theta2, theta3)

    x = cos(theta1) + cos(theta1 + theta2) + cos(theta1 + theta2 + theta3);
    y = sin(theta1) + sin(theta1 + theta2) + sin(theta1 + theta2 + theta3);
    phi = theta1 + theta2 + theta3;
end

