%% Channel DNS Subfunction - Correct
%% Purpose
%   Calculate the pressure and velocity correct field
%% Method
%   The correct pressure field satisfies:
%
%       laplace P_c+ = 0
%       P_c+(y=1) = 1, P_c+(y=-1) = 0                   (Selected)
%
%       laplace P_c- = 0
%       P_c-(y=1) = 0, P_c-(y=-1) = 1
%
%   And the correct velocity field satisfies:
%
%       laplace V_c+ - V_c+ / niu*dt = 1/niu * P_c+
%       V_c+(y = -1 and 1) = 0                          (Selected)
%
%       laplace V_c- - V_c- / niu*dt = 1/niu * P_c-
%       V_c-(y = -1 and 1) = 0
%
%   Considering the symmetry of P_c+ and P_c-, V_c+ and V_c-,
%   we only need to compute one pair of the equations. In this
%   function, the first pair of equations are solved.
%% Parameters
%   Input parameters:
%   nx ----------------------- node number in x
%   ny ----------------------- node number in y
%   nz ----------------------- node number in z
%   alpha -------------------- wave number in x
%   beta --------------------- wave number in z
%   re ----------------------- reynold number
%   dt ----------------------- time space
%   Output parameters:
%   pc ----------------------- correct pressure field
%   uc ----------------------- correct u-velocity field
%   vc ----------------------- correct v-velocity field
%   wc ----------------------- correct w-velocity field
%   M ------------------------ correct matrix in physical
%   Note : All output parameters are in spectral space
%% Author
%   Written by Luohao Wang on 2015-9-10
%   Contact : lh-wang13@mails.tsinghua.edu.cn

%% Code
function [pc, uc, vc, wc, M] = correct(nx, ny, nz, alpha, beta, re, dt)
pc = zeros(nx, ny, nz);
uc = zeros(nx, ny, nz);
vc = zeros(nx, ny, nz);
wc = zeros(nx, ny, nz);
M = zeros(nx, nz, 2, 2);
s  = zeros(1, ny);
gamma0 = 11/6;
%Pressure correct
for iter_x = 0:nx-1
    for iter_z = 0:nz-1
        k = alpha^2*(iter_x - nx/2)^2 + beta^2*(iter_z - nz/2)^2;
        pc(iter_x+1,:,iter_z+1) = ode2_solve(k, s, 1, [0,1], 'p', 's');
    end
end
%Velocity correct
d_p_dx = globe_partial(pc, 'x', alpha);
d_p_dy = globe_partial(pc, 'y', 0);
d_p_dz = globe_partial(pc, 'z', beta);
for iter_x = 0:nx-1
    for iter_z = 0:nz-1
        uc(iter_x+1,:,iter_z+1) = ode2_solve(gamma0/dt, d_p_dx(iter_x+1,:,iter_z+1), 1, [0,0], 's', 's');
        vc(iter_x+1,:,iter_z+1) = ode2_solve(gamma0/dt, d_p_dy(iter_x+1,:,iter_z+1), 1, [0,0], 's', 's');
        wc(iter_x+1,:,iter_z+1) = ode2_solve(gamma0/dt, d_p_dz(iter_x+1,:,iter_z+1), 1, [0,0], 's', 's');
    end
end
uc = re * uc;
vc = re * vc;
wc = re * wc;
%form correct matrix
div_V = globe_partial(uc, 'x', alpha) + globe_partial(vc, 'y', 0) ...
    + globe_partial(wc, 'z', beta);
div_V = y_trans(div_V, -1);
for iter_x = 1:nx
    for iter_z = 1:nz
        M(iter_x, iter_z, 1, 1) = div_V(iter_x, 1, iter_z);
        M(iter_x, iter_z, 1, 2) = div_V(iter_x, ny, iter_z);
        M(iter_x, iter_z, 2, 1) = div_V(iter_x, ny, iter_z);
        M(iter_x, iter_z, 2, 2) = div_V(iter_x, 1, iter_z);
    end
end
%transform to physical space
pc = y_trans(pc, -1);
uc = y_trans(uc, -1);
vc = y_trans(vc, -1);
wc = y_trans(wc, -1);
end