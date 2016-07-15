%test fcl
u = ones(32,32,32);
for x = 1:32
    for z = 1:32
        u(x,1,z) = ((x-16)/32)^2;%sin((x-1)/32*2*pi);
    end
end
u(:,2:32,:) = 0;
% v = u;
% for iter = 1:32
%     temp = reshape(u(:,iter,:), 32,32);
%     temp = fft(temp);
%     temp = fftshift(temp,1);
%     temp = fft(temp.');
%     temp = fftshift(temp,1);
%     us(:,iter,:) = temp.';
%     %       us(:,iter,:) = fftshift(fft(u(:,iter,:)));
% end
% for x = 1:32
%     for z = 1:32
%         us(x,:,z) = FCT(us(x,:,z),1);
%     end
% end
us = global_trans(u,1);
%  dus = globe_partial(us, 'x',1);
%  dus = globe_partial(dus, 'x',1);
% dus = globe_partial(dus, 'x',1);
% us = globe_partial(dus, 'x',1);
% vs = us;
% ws = fcl(us,vs);
% for x = 1:32
%     for z = 1:32
%         ws(x,:,z) = FCT(ws(x,:,z),-1);
%     end
% end
%  w = zeros(size(u));
% for iter = 1:32
%     temp = reshape(ws(:,iter,:),32,32);
%     temp = ifftshift(temp,1);
%     temp = ifft(temp);
%     temp = temp.';
%     temp = ifftshift(temp,1);
%     temp = ifft(temp).';
%     w(:,iter,:) = temp;
% end
w = global_trans(us,-1);
err = w - u;
ww = real(reshape(w(:,1,:),32,32));
surf(ww);
%err = ww - reshape(u(:,1,:).^2,32,32);
%surf(err);