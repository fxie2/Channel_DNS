% %test slice
% n = 100;
% a = zeros(n,n,n);
% for x = 1:n
%     for y = 1:n
%         for z = 1:n
%             a(x,y,z) = exp(-((x/100 - 0.5)^2 + (y/100-0.5)^2+(z/100-0.5)^2));
%         end
%     end
% end
% slice(a,n/2,n/2,n/2);
% n=3;   %离散点数
% 
% %随机生成点坐标和矢量值
% x=rand(n,1);
% y=rand(n,1);
% z=rand(n,1);
% u=rand(n,1);
% v=rand(n,1);
% w=rand(n,1);
% 
% quiver3(x,y,z,u,v,w)     %画矢量
% hold on
% plot3(x,y,z,'o'), hold off         %画点
% 
% xlabel('X Axis')
% ylabel('Y Axis')
% zlabel('Z Axis')
