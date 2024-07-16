% Define the error function
function error = error_function(x_desired, y_desired, phi_desired, theta1, theta2, theta3)
    % This function calculates the error between the desired end-effector
    % position (x_desired, y_desired, phi_desired) and the actual
    % end-effector position calculated using forward kinematics based on
    % the joint angles (theta1, theta2, theta3) 
    
    [x, y, phi] = forward_kinematics(theta1, theta2, theta3);
    error = sqrt((x_desired - x)^2 + (y_desired - y)^2 + (phi_desired - phi)^2);
end