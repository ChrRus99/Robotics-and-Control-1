% Define the Jacobian function
function J = jacobian(theta1, theta2, theta3)
    % This function return the Jacobian matrix J at the given joint angles
    % (theta1, theta2, theta3) 
    
    J = [-sin(theta1)-sin(theta1+theta2)-sin(theta1+theta2+theta3), -sin(theta1+theta2)-sin(theta1+theta2+theta3), -sin(theta1+theta2+theta3); 
         cos(theta1)+cos(theta1+theta2)+cos(theta1+theta2+theta3), cos(theta1+theta2)+cos(theta1+theta2+theta3), cos(theta1+theta2+theta3)
         1, 1, 1];
end
