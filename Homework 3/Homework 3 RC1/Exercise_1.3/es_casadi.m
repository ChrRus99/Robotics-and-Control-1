import casadi.*

% define a 3x3 matrix
A = SX.sym('a', [3,3])

% define 3x1 vector
b = SX.sym('b', [3,1])

% compute the matrix vector product c = A*b
c = A*b

% apply sin funcion to the first element of d = A*b
d = A*b;
d(1) = sin(d(1))

% compute derivatives
diff_c_b = jacobian(c,b)
diff_d1_b = gradient(d(1),b)

% transform expressions in functions
f_c = Function('f_c',{A,b},{c});
f_diff_d1_b = Function('f_diff_d1_b',{A,b},{diff_d1_b});

% evaluate the functions
f_c(eye(3), [1,2,3])
f_diff_d1_b(eye(3), [1,2,3])

% generate the mex functions
opts = struct('main', true,...
              'mex', true);
f_c.generate('f_c_mex.c',opts);
mex f_c_mex.c
f_diff_d1_b.generate('f_diff_d1_b_mex.c',opts);
mex f_diff_d1_b_mex.c

% call the mex function
f_c_mex('f_c',eye(3), [1,2,3])
f_diff_d1_b_mex('f_diff_d1_b',eye(3), [1,2,3])