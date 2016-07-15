%% Channel DNS Subfunction - partial
%% Purpose
%   Calculate the partial derivative in spectral space
%% Parameters
%   Input parameters:
%   f ------------------------- input data in spectral space
%   direction ----------------- 'x' for d_dx, and so is others
%   parameter ----------------- alpha for 'x', beta for 'z', 
%   nothing for 'y'
%% Author
%   Written by Luohao Wang on 2015-9-11
%   Contact : lh-wang13@mails.tsinghua.edu.cn

%% Code
function d = partial(f, direction, parameter)
len = length(f);
if direction == 'x'     %now parameter == alpha
    d = complex(0,1)*parameter.*(-len/2:len/2-1).*f;
elseif direction == 'z' %now parameter == beta
    d = complex(0,1)*parameter.*(-len/2:len/2-1).*f;
elseif direction == 'y' %now parameter is useless, can be anything
    P = len - 1;
    [~,N] = size(f);
    if N == 1
        f = transpose(f);
    end
    d = zeros(size(f));
    d(1) = sum((1:2:P).*f((1:2:P)+1));
    for p = 2:P
        d(p) = 2*sum((p:2:P).*f((p:2:P)+1));
    end
    if N == 1
        d = transpose(d);
    end
end
end