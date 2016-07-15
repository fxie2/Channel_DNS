%% Channel DNS Subfunction - global_trans
%% Purpose
%   Transform between physical space and spectral space
%% Parameters
%   Input parameters:
%   u ------------------------ spectral or physical data
%   direction ---------------- 1 for p -> s, -1 for s -> p
%   Output parameter:
%   w ------------------------ physical or spectral result
%% Author
%   Written by Luohao Wang on 2015-9-13
%   Contact : lh-wang13@mails.tsinghua.edu.cn

%% Code
function w = global_trans(u, direction)
[nx, ny, nz] = size(u);
w = zeros(size(u));
if direction == 1
    for iter_y = 1:ny
        temp = reshape(u(:,iter_y,:), nx, nz);
%         temp = fft(temp);
%         temp = fftshift(temp, 1);
%         temp = fft(temp.');
%         temp = fftshift(temp, 1);
%         w(:,iter_y,:) = temp.';
        temp = fftshift(fft2(temp));
        w(:,iter_y,:) = temp;
    end
    for iter_x = 1:nx
        for iter_z = 1:nz
            w(iter_x,:,iter_z) = FCT(w(iter_x,:,iter_z), 1);
        end
    end
elseif direction == -1
    for iter_x = 1:nx
        for iter_z = 1:nz
            w(iter_x,:,iter_z) = FCT(u(iter_x,:,iter_z), -1);
        end
    end
    for iter_y = 1:ny
        temp = reshape(w(:,iter_y,:), nx, nz);
%         temp = ifftshift(temp, 1);
%         temp = ifft(temp).';
%         temp = ifftshift(temp,1);
%         temp = ifft(temp).';
%         w(:,iter_y,:) = temp;
        temp = ifft2(ifftshift(temp));
        w(:,iter_y,:) = temp;
    end
    w = real(w);
end
end