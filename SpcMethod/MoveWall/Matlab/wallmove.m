%% Channel DNS Subfunction - wallmove
%% Purpose
%   Generate physical shape of both upper and lower wall
%   with form:
%
%   etau(x,z,t) = sum Ax_n*sin(kx_n*x - cx_n*t + phi0x_n)
%                 sum Az_n*sin(kz_n*z - cz_n*t + phi0z_n), n = 1,2,...
%
%   etad(x,z,t) has the same form.
%% Parameters
%   Input parameters:
%   nx, nz ---------------------------- node number in x, z
%   alpha, beta ----------------------- wave number in x, z
%   Ax, Az ---------------------------- wave amptitude
%   kx, kz ---------------------------- wall wave number
%   cx, cz ---------------------------- wall wave phase spead
%   phi0x, phi0z ---------------------- wall wave initial phase
%   Output parameters:
%   etau, etad ------------------------ physical wall movement
%% Author
%   Written by Luohao Wang on 2015-9-28
%   Contact : lh-wang13@mails.tsinghua.edu.cn

%% Code
function [etau, etad] = wallmove(nx, nz, alpha, beta, Ax, kx, cx, phi0x, ...
    Az, kz, cz, phi0z, t)
x = linspace(0, 2*pi/alpha, nx+1);
x = x(1:nx);
z = linspace(0, 2*pi/beta, nz+1);
z = z(1:nz);
etau = zeros(nx, nz);
etad = zeros(nx, nz);
[len1, len2] = size(Ax);
[len3, len4] = size(Az);
if (len1 ~= 4 || len3 ~= 4)
    disp('Dimension error for A in wallmove function!');
else
    %calculate in x direction
    for iter_z = 1:nz
        for iter_wave = 1:len2
            etau(:,iter_z) = etau(:,iter_z) + Ax(iter_wave, 1)*...
                sin(kx(iter_wave, 1)*x - cx(iter_wave, 1)*t + phi0x(iter_wave, 1));
            etad(:,iter_z) = etad(:,iter_z) + Ax(iter_wave, 2)*...
                sin(kx(iter_wave, 2)*x - cx(iter_wave, 2)*t + phi0x(iter_wave, 2));
        end
    end
    %calculate in z direction
    for iter_x = 1:nx
        for iter_wave = 1:len4
            etau(iter_x,:) = etau(iter_x,:) + Az(iter_wave, 1)*...
                sin(kz(iter_wave, 1)*z - cz(iter_wave, 1)*t + phi0z(iter_wave, 1));
            etad(iter_x,:) = etad(iter_x,:) + Az(iter_wave, 2)*...
                sin(kz(iter_wave, 2)*z - cz(iter_wave, 1)*t + phi0z(iter_wave, 1));
        end
    end
end
end