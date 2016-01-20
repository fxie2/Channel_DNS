%% Channel DNS Subfunction - globe_partial
%% Purpose
%   Calculate the parital derivative of the hole field
%% Parameters
%   Input parameters:
%   f ------------------------- spectral data in 3-D
%   direction ----------------- 'x' for d_dx, and so on
%   parameter ----------------- alpha for 'x', beta for 'z', nothing for
%   'y'
%% Author
%   Written by Luohao Wang on 2015-9-12
%   Contact : lh-wang13@mails.tsinghua.edu.cn

%% Code
function df = globe_partial(f, direction, parameter)
[nx, ~, nz] = size(f);
df = zeros(size(f));
if direction == 'x'
    for iter_x = 0:nx-1
        df(iter_x+1,:,:) = complex(0,1).*parameter.*f(iter_x+1,:,:).*(iter_x - nx/2);
    end
elseif direction == 'z'
    for iter_z = 0:nz-1
        df(:,:,iter_z+1) = complex(0,1).*parameter.*f(:,:,iter_z+1).*(iter_z - nz/2);
    end
elseif direction == 'y'
    for iter_x = 1:nx
        for iter_z = 1:nz
            df(iter_x,:,iter_z) = partial(f(iter_x,:,iter_z), 'y', 0);
        end
    end
else
    disp('direction error');
end
end