% %test directive
% x = linspace(0,2*pi,1001);
% x = x(1:1000)';
% u = sin(x) + cos(x);
% us = fft(u);
% us = fftshift(us);
% dus = us.*complex(0,1).*(-500:499)';
% du = ifftshift(dus);
% du = ifft(du);
% plot(x, du);
u = fftshift(fft(sin([0:255]/256*2*pi)));
v = u;
w = fcl(u,v);
a = ifft(ifftshift(w));
err = a - sin([0:255]/256*2*pi).^2;
plot([0:255],err)