%% Channel DNS Subroutine - dpdx_correct
%% Puporse
%   Correct the gradient of pressure according to the constant mass flow
%   boundary condition
%% Method
%   Treat dpdx as the function of mass flow rate, use secant method to find
%   the sloution

%% Code
function dpdx = dpdx_correct(u, dpdx_u, u_hist, dpdx_hist, q_set, beta)
[nx, ny, nz] = size(u);
u_p = global_trans(u, -1);
u_hist_p = global_trans(u_hist, -1);
u_diff = u_p(nx, :, :) - u_hist_p(nx, :, :);
u_diff = reshape(u_diff, ny, nz);
if ~(max(max(abs(u_diff))) >= 1e-6 && abs(dpdx_u-dpdx_hist) > eps())
    dpdx = dpdx_u;
else
    dz = 2*pi/beta/nz;
    q_u = 0;
    q_hist = 0;
    for iter_z = 1:nz
        u_ys_temp = reshape(FCT(u_p(nx, :, iter_z), 1), 1, ny);
        u_hist_ys_temp = reshape(FCT(u_hist_p(nx, :, iter_z), 1), 1, ny);
        for iter_y = 1:2:ny
            n = iter_y - 1;
            q_u = q_u + u_ys_temp(iter_y) * 2 / (1-n^2) * dz;
            q_hist = q_hist + u_hist_ys_temp(iter_y) * 2 / (1-n^2) * dz;
        end
    end
    if abs(q_u - q_hist) <= 1e-6
        dpdx = dpdx_u;
    else
        dpdx = dpdx_u + (q_set - q_u) * (dpdx_u - dpdx_hist) / (q_u - q_hist);
    end
end
end