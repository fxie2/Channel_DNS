%% Channel DNS Subfunction - initialize
%% Purpose
%   Initialize the fluid field at the start time
%% Parameters
%   nx, ny, nz ----------------------------- mesh size
%   alpha, beta ---------------------------- wave number in x and z
%% Author
%   Written by Luohao Wang on 2015-9-12
%   Contact : lh-wang13@mails.tsinghua.edu.cn

%% Code
function [p, u, v, w] = initialize(nx, ny, nz, alpha, dpdx)
p = zeros(nx, ny, nz);
u = zeros(nx, ny, nz);
v = zeros(nx, ny, nz);
w = zeros(nx, ny, nz);

%Use theoretical Poiseuille flow field to initialize

%Velocity initialize
for iter_y = 1:ny
    y = cos((iter_y-1)*pi/(ny-1));
    u(:, iter_y, :) = 1 - y^2;
end
%Pressure initialize
for iter_x = 2:nx
    p(iter_x,:,:) = p(iter_x-1,:,:) - dpdx*2*pi/alpha/nx;
end
end