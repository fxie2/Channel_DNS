%% Channel DNS Subfunction - Gterm
%% Purpose
%   Calculate the G term in Pressure step
%% Parameters
%   Input parameters:
%   fy -------------------------- lamb-y term
%   us, vs, ws ------------------ velocity term
%   alpha, beta ----------------- wave number
%   re -------------------------- reynold number
%   Output parameter
%   G --------------------------- G term
%% Author
%   Written by Luohao Wang on 2015-9-12
%   Contact : lh-wang13@mails.tsinghua.edu.cn

%% Code
function G = Gterm(fy, us, vs, ws, alpha, beta, re)
rot_x = globe_partial(ws, 'y', 0) - globe_partial(vs, 'z', beta);
rot_z = globe_partial(vs, 'x', alpha) - globe_partial(us, 'y', 0);
G = fy;
[nx, ~, nz] = size(fy);
i = complex(0,1);
parfor iter_x = 1:nx
    for iter_z = 1:nz
        G(iter_x,:,iter_z) = G(iter_x,:,iter_z) + (i*alpha*(iter_x-1-nx/2)*rot_z(iter_x,:,iter_z)...
            - i*beta*(iter_z-1-nz/2)*rot_x(iter_x,:,iter_z)) / re;
    end
end
end