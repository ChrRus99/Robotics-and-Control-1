%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SCARA inverse kinematics via inverse differential kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% import functions and define the target trajectory

% % import casadi
% import casadi.*
% % import kin and J functions (to use casadi function instead of mex)
% load f_J.mat
% load f_x.mat

% initial configuration
q_0 = [-pi/4,pi/4,-0.15,0];
x_0 = f_x_mex('f_x',q_0);

% define time
T_sampling = 0.0001;
T_horizon = 2;
omega = 2*2*pi/T_horizon; 
t = 0:T_sampling:T_horizon;
num_t = length(t);

% define the end-effector trajectory
r = 0.2;
c_x = x_0(1)-r;
c_y = x_0(2);
c_z = x_0(3);
a_z = 0.3/2;
x_target = [c_x + r*cos(omega*t);
            c_y + r*sin(omega*t);
            c_z + a_z*t;
            zeros(1,num_t)]';
x_dot_target = [-omega*r*sin(omega*t);
                 omega*r*cos(omega*t);
                 a_z*ones(1,num_t);
                 zeros(1,num_t)]';

% plot target trajectory
figure()
plot3(x_target(:,1), x_target(:,2), x_target(:,3))
hold on
plot3(x_target(1,1), x_target(1,2), x_target(1,3), 'X', 'color', 'r', 'LineWidth', 2)
plot3(x_target(end,1), x_target(end,2), x_target(end,3), 'X', 'color', 'g', 'LineWidth', 2)
hold on
legend('P x target', 'starting point', 'ending point')
grid on
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')


%% compute the joint trajectory

% allocate the space
q_d = zeros(num_t,4);
q_d_dot = zeros(num_t,4);

% set then first joint configuration
q_d(1,:) = q_0;

% compute q_target_euler
for i = 2 : num_t
    % compute J (with mex)
    J = f_J_mex('f_J', q_d(i-1,:));
%     % compute J (with casadi function)
%     J = full(f_J(q_d(i-1,:))); 
    % compute q_dot_target
    q_d_dot(i-1,:) = (eye(4)/J)*(x_dot_target(i-1,:)');
    % compute q_target
    q_d(i,:) = q_d(i-1,:) + T_sampling*q_d_dot(i-1,:);
end



%% Compute the end-effector trajectory performed

% allocate the space
x_inv_diff = zeros(num_t,4);

% get x_e using f_x
for i = 1 : num_t
   % using mex function 
   x_inv_diff(i,:) = f_x_mex('f_x', q_d(i,:));
%    % using casadi function
%    x_inv_diff(i,:) = full(f_x(q_d(i,:))); 
end

% compute the error
RMSE = sqrt(mean(mean((x_target-x_inv_diff).^2)))

% compare target trajectory with inv diff kin trj
figure()
plot3(x_target(:,1), x_target(:,2), x_target(:,3))
hold on
plot3(x_inv_diff(:,1), x_inv_diff(:,2), x_inv_diff(:,3))
legend('P x target', 'P x inv diff kin')
grid on
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')