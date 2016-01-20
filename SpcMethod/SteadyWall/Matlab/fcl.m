%% Channel DNS Subfunction - fcl
%% Purpose
%   Compute u_i*v_i in spectral space using FFT instead of convolution
%% Method
%   See 3/2 rule
%% Parameters
%   Input parameters:
%   u, v ------------------------ spectral data in 3-D for convolution
%   Output parameter:
%   w --------------------------- spectral result in 3-D
%% Attention
%   For 3-D reason, we cannot directly perform fft2 or ifft2 on the slice
%   like u(:,iter,:) which will lead to dimension error. So we have to use
%   a temporary matrix to store the slice, mannuly do twice fft and
%   fftshift, and return the result to the initial slice.
%% Author
%   Written by Luohao Wang on 2015-9-12
%   Contact : lh-wang13@mails.tsinghua.edu.cn

%% Code
function w = fcl(u, v)
[nx, ny, nz] = size(u);
U_trans_y = zeros(size(u));
V_trans_y = zeros(size(v));
for iter_x = 0:nx-1
    for iter_z = 0:nz-1
        U_trans_y(iter_x+1,:,iter_z+1) = FCT(u(iter_x+1,:,iter_z+1),-1);
        V_trans_y(iter_x+1,:,iter_z+1) = FCT(v(iter_x+1,:,iter_z+1),-1);
    end
end
U_extend = zeros(nx/2*3, ny, nz/2*3);
V_extend = zeros(nx/2*3, ny, nz/2*3);
U_extend(nx/4+1:nx/4+nx,:,nz/4+1:nz/4+nz) = U_trans_y;
V_extend(nx/4+1:nx/4+nx,:,nz/4+1:nz/4+nz) = V_trans_y;
U_real = zeros(size(U_extend));
V_real = zeros(size(V_extend));
for iter_y = 0:ny-1
    temp = reshape(U_extend(:,iter_y+1,:), nx/2*3, nx/2*3);
%     temp = ifftshift(temp,1);
%     temp = ifft(temp).';
%     temp = ifftshift(temp,1);
%     temp = ifft(temp).';
    temp = ifft(ifftshift(temp));
    U_real(:,iter_y+1,:) = temp;
    temp = reshape(V_extend(:,iter_y+1,:), nx/2*3, nx/2*3);
%     temp = ifftshift(temp,1);
%     temp = ifft(temp).';
%     temp = ifftshift(temp,1);
%     temp = ifft(temp).';
    temp = ifft(ifftshift(temp));
    V_real(:,iter_y+1,:) = temp;   %temp here is only for coefficient shift
end
W_real = U_real.*V_real;
W_extend = zeros(size(U_extend));
for iter_y = 0:ny-1
    temp = reshape(W_real(:,iter_y+1,:), nx/2*3, nx/2*3);
%     temp = fft(temp);
%     temp = fftshift(temp,1).';
%     temp = fft(temp);
%     temp = fftshift(temp,1).';
    temp = fftshift(fft(temp));
    W_extend(:,iter_y+1,:) = temp;
end
W_trans_y = W_extend(nx/4+1:nx/4+nx,:,nz/4+1:nz/4+nz);
w = zeros(size(u));
for iter_x = 0:nx-1
    for iter_z = 0:nz-1
        w(iter_x+1,:,iter_z+1) = FCT(W_trans_y(iter_x+1,:,iter_z+1),1);
    end
end
w = 9/4*w;
end