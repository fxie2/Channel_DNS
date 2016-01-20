%% Channel DNS Subfunction - lamb
%% Purpose
%   Calculate the lamb term F = V X w.
%% Parameters
%   Input parameters:
%   us , vs, ws ------------------------- spectral velocity in 3-D
%   alpha, beta ------------------------- wave number in x and z
%   Output parameters:
%   fx, fy, fz -------------------------- convective term in spectral
%% Author
%   Written by Luohao Wang on 2015-9-12
%   Contact : lh-wang13@mails.tsinghua.edu.cn

%% Code
function [fx, fy, fz] = lamb(us, vs, ws, alpha, beta)
rot_x = globe_partial(ws, 'y', 0) - globe_partial(vs, 'z', beta);
rot_y = globe_partial(us, 'z', beta) - globe_partial(ws, 'x', alpha);
rot_z = globe_partial(vs, 'x', alpha) - globe_partial(us, 'y', 0);
fx = fcl(vs, rot_z) - fcl(ws, rot_y);
fy = fcl(ws, rot_x) - fcl(us, rot_z);
fz = fcl(us, rot_y) - fcl(vs, rot_x);
end