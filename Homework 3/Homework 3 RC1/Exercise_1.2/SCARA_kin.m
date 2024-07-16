%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Script that computes the SCARA forward kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set/import variables

% import casadi
import casadi.*

% define the symbolic variables
num_dof = 4;
q = SX.sym('q', [num_dof,1]);

% define the DH parameters
alpha_list = [0, 0, 0, 0];
a_list = [0.5, 0.4, 0, 0.2];
d_list = [0.5, -0.005, -0.15, -0.005];
theta_list = [0, 0, 0, 0];

% define a vector with the joint type (if true the joint is revolute)
joint_type_list = [true, true, false, true];



%% compute the DH relative transformations

% allocate the space
R_DH = cell(num_dof,1);
l_DH = cell(num_dof,1);

%iterate all the joints
for j = 1 : num_dof
    
    % check the joint type and determine d_j and theta_j
    if joint_type_list(j)
        % if the joint is rev theta_j = theta_j^0 + q_j
        theta_j = theta_list(j) + q(j);
        d_j = d_list(j);
    else
        % if the joint is prismatic d_j = d_j^0 + q_j
        theta_j = theta_list(j);
        d_j = d_list(j) + q(j);
    end
    
    % get DH transformation
    [R_DH{j}, l_DH{j}] = DH_transform(alpha_list(j), a_list(j), d_j, theta_j);
    
end



%% get orientation and translation of the links w.r.t. world

% allocate variables
R_j_0_list = cell(num_dof+1,1);
l_j_0_list = cell(num_dof+1,1);

% initialize the variables
l_j_0_list{1} = SX.zeros(3,1);
R_j_0_list{1} = SX.eye(3);

%compute the relative transformations
for j = 2 : num_dof+1
   % l_j_0 = l_j-1_0 + R_j-1^0*l_j^j-1
   l_j_0_list{j} = l_j_0_list{j-1} + R_j_0_list{j-1}*l_DH{j-1};
   % R_j_0 = R_j-1^0*R_j^j-1
   R_j_0_list{j} = R_j_0_list{j-1}*R_DH{j-1};
end


%% convert R_e_0 in yaw pitch and roll angles

[yaw, pitch, roll] = rot_2_YPR_angles(R_j_0_list{num_dof+1});



%% get the numeric function and the mex function

% get functions evaluable numerically
l_matrix = [l_j_0_list{2},...
            l_j_0_list{3},...
            l_j_0_list{4},...
            l_j_0_list{5}];
R_matrix = [R_j_0_list{2},...
            R_j_0_list{3},...
            R_j_0_list{4},...
            R_j_0_list{5}];
f_kin_l = Function('f_kin_l',{q}, {l_matrix});
f_kin_R = Function('f_kin_R',{q}, {R_matrix});
f_x = Function('f_x',{q},{[l_j_0_list{num_dof+1}; yaw]});

% save the casadi function
save('f_kin_l','f_kin_l')
save('f_kin_R','f_kin_R')
save('f_x','f_x')

% generate the mex function
opts = struct('main', true,...
              'mex', true);
f_kin_l.generate('f_kin_l_mex.c',opts);
f_kin_R.generate('f_kin_R_mex.c',opts);
f_x.generate('f_x_mex.c',opts);
mex f_kin_l_mex.c
mex f_kin_R_mex.c
mex f_x_mex.c


%% test the kinematic

% set a configuration
q_num = [0,0,0,0];

% % get the forward kin numeric values (usign casadi function)
% l_matrix_num = full(f_kin_l(q_num));
% R_matrix_num = full(f_kin_R(q_num));
% x = full(f_x(q_num));

% get the forward kin numeric values (usign mex function)
l_matrix_num = f_kin_l_mex('f_kin_l', q_num);
R_matrix_num = f_kin_R_mex('f_kin_R', q_num);
x = f_x_mex('f_x', q_num);

% visualize a schematic representation of the robot
print_SCARA_kin(l_matrix_num, R_matrix_num)



%% derive the Jacobian

% compute the jacobian
J = jacobian(f_x(q), q);

% get the casadi function
f_J = Function('f_J',{q}, {J});

% save the casadi function
save('f_J','f_J')

% generate the mex function
f_J.generate('f_J_mex.c',opts);
mex -largeArrayDims f_J_mex.c



%% get the jacobian in q = [0,0,0,0]

% using mex function
J_0 = full(f_J_mex('f_J',[0,0,0,0]))
% % using casadi function
% J_0 = full(f_J([0,0,0,0]))