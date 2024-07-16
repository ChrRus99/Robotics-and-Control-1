%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SCARA Derivation of the dynamics equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% import casadi and define robot parameters
import casadi.*

% define the DH parameters
alpha_list = [0, 0, 0, 0];
a_list = [0.5, 0.4, 0, 0.2];
d_list = [0.5, -0.005, -0.15, -0.005];
theta_list = [0, 0, 0, 0];

% links mass
m_vect = [25,20,15,10,5];

% links length
l_vect = [0.5,0.5,0.4,0.3, 0.2];
r_link = 0.1;

% define the links center of mass coordinates
% (expressed in the correspondent DH reference frame)
p_matrix = [ 0.,  0., l_vect(1)/2;
            -l_vect(2)/2, 0., 0.;
            -l_vect(3)/2, 0., 0.;
             0.,  0., l_vect(4)/2;
             -l_vect(5)/2, 0., 0.];
         
% define the links zz inertia
% (expressed in the correspondent DH reference frame)
I_zz_vect = [m_vect(1)/12*(l_vect(1)^2+3*r_link^2),...
             m_vect(2)/12*(l_vect(2)^2+3*r_link^2),...
             m_vect(3)/12*(l_vect(3)^2+3*r_link^2),...
             0.5*m_vect(4)*r_link^2,...
             m_vect(5)/12*(l_vect(5)^2+3*r_link^2)];

% define a vector with the joint type (if true the joint is revolute)
num_dof = 4;
joint_type_list = [true, true, false, true];



%% define the symbolic variables

% joint positions
q = SX.sym('q', [num_dof,1]);
% joint coordinates
q_dot = SX.sym('q_dot', [num_dof,1]);



%% compute the expressions of the DH reference frames

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



%% compute the positions of the centers of mass

% allocate variables
P_p = cell(num_dof+1,1);

% compute positions
for j = 1 : num_dof+1
    P_p{j} = l_j_0_list{j} + R_j_0_list{j}*p_matrix(j,:)';
end



%% express Z-axis orientation using the yaw angle

% allocate variables
yaw_angles = cell(num_dof+1,1);
alpha_vect{1} = 0;

% compute Z angle
for j = 2 : num_dof+1
    [yaw_angles{j}, pitch, roll] = rot_2_YPR_angles(R_j_0_list{j});
end



%% get jacobians

% allocate the space
J_p = cell(num_dof+1,1);
J_o = cell(num_dof+1,1);

for j = 2 : num_dof+1
    % linear velocities
    J_p{j} = jacobian(P_p{j}, q);
    % angular velocties
    J_o{j} = jacobian(yaw_angles{j}, q);
end



%% compute the inertia matrix and the potential energy

% init the matrix with the first link
B = m_vect(2)*J_p{2}'*J_p{2} + I_zz_vect(2)*J_o{2}'*J_o{2};
g = [0;0;-9.80665];
U = -m_vect(2)*g'*P_p{2};

for j = 3 : num_dof+1
    B = B + m_vect(j)*J_p{j}'*J_p{j} + I_zz_vect(j)*J_o{j}'*J_o{j};
    U = U -m_vect(j)*g'*P_p{j};
end 

%% compute dynamics equations

% terms dependent on velocities
H = cell(num_dof,1);
for i=1:num_dof
    % get H1 and H2
    H{i} = jacobian(B(i,:),q) - 0.5*reshape(jacobian(B,q(i)), num_dof, num_dof);
end
C_vect = [q_dot'*H{1}*q_dot;
          q_dot'*H{2}*q_dot;
          q_dot'*H{3}*q_dot;
          q_dot'*H{4}*q_dot];

% gravitational contribution
g_vect = jacobian(U, q);



%% get functions evaluable numerically

% get casadi functions
f_B = Function('f_B',{q}, {B});
f_U = Function('f_U',{q}, {U});
f_C_vect = Function('f_C_vect',{q, q_dot}, {C_vect});
f_g_vect = Function('f_g_vect',{q}, {g_vect});

% get mex functions
opts = struct('main', true,...
              'mex', true);
f_B.generate('f_B_mex.c',opts);
mex f_B_mex.c
f_U.generate('f_U_mex.c',opts);
mex f_U_mex.c
f_C_vect.generate('f_C_vect_mex.c',opts);
mex f_C_vect_mex.c
f_g_vect.generate('f_g_vect_mex.c',opts);
mex f_g_vect_mex.c



%% Test gravity

% expected gravitational contribution
F_g_expected = (m_vect(4)+m_vect(5))*g

% gravitational contribution computed by the model
F_g_0 = -full(f_g_vect_mex('f_g_vect', [0,0,0,0]))
F_g_1 = -full(f_g_vect_mex('f_g_vect', [0,0,0.1,0]))



%% Test inertia matrix

% compute the inertia matrix in [0,0,0,0]
B_0 = full(f_B_mex('f_B',[0,0,0,0]))

% expected value of B_33
B_33_expected = m_vect(4)+m_vect(5)

% expected value of B_44
B_44_expected = p_matrix(5,1)^2*m_vect(5) + I_zz_vect(5)

% compute the inertia matrix in [0,0,0,0]
B_1 = full(f_B_mex('f_B',[pi/4,-pi/4,0,0]))



%% Test fictitious contributions

% compute fictitious contributions with
% q=[0,pi/2,0,0] and q_dot=[1,0,0,0]
C_vect_1 = full(f_C_vect_mex('f_C_vect', [0,pi/2,0,0], [1,0,0,0]))

% compute fictitious contributions with
% q=[0,pi/2,0,0] and q_dot=[2,0,0,0]
C_vect_2 = full(f_C_vect_mex('f_C_vect', [0,pi/2,0,0], [2,0,0,0]))

% compute fictitious contributions with
% q=[0,pi/2,0,0] and q_dot=[3,0,0,0]
C_vect_3 = full(f_C_vect_mex('f_C_vect', [0,pi/2,0,0], [3,0,0,0]))