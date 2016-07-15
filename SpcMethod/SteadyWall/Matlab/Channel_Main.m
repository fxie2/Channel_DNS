%% Channel Flow DNS Program
%% Purpose
%   Simulate the fluid field of 3-D channel flow
%% Method
%   DNS with spectral method. Fouries transform in x and z direction,
%   Chebyshv transform in y direction
%   Galerkin-Tau Method
%% Author
%   Written by Luohao Wang on 2015-9-10
%   Contact : lh-wang13@mails.tsinghua.edu.cn

%% Basic parameters
nx = 8;   %node number in x
ny = 17;   %node number in y
nz = 8;   %node number in z
nxh = nx/2; %half of nx
nyh = ny/2; %half of ny
nzh = nz/2; %half of nz
nx32 = nx *3/2; %3/2 rule
ny32 = ny *3/2; %3/2 rule
nz32 = nz *3/2; %3/2 rule
re = 1;       %reynold number
alpha = 1.0;    %wave number in x
beta = 1.0;     %wave number in z
dt = 0.01;         %time space
start_time = 0; %start time
end_time = dt*100;   %end time
step_num = round((end_time - start_time) / dt); %time steps

%% Boundary conditions
dpdx = 2/re;                %Constant pressure gradient boundary condition
mass_flux = 8/3*pi/beta;              %Constant mass flux condition, 0 for disabled

%% Fluid Field Data
% x = linspace(0,2*pi/alpha,nx+1);  %x coordinate
% x = x(1:nx);
y = cos((0:(ny-1))/(ny-1)*pi);            %y coordinate
% z = linspace(0,2*pi/beta,nz+1);   %z coordinate
% z = z(1:nz);
% p = zeros(nx, ny, nz);          %physical total pressure
% u = zeros(nx, ny, nz);          %physical u-velocity
% v = zeros(nx, ny, nz);          %physical v-velocity
% w = zeros(nx, ny, nz);          %physical w-velocity
[p, u, v, w] = initialize(nx, ny, nz, alpha, dpdx);
% v = 1e-14*randn(size(v));
u0 = u;
ps = global_trans(p, 1);        %spectral total pressure
us = global_trans(u, 1);        %spectral u-velocity
vs = global_trans(v, 1);        %spectral v-velocity
ws = global_trans(w, 1);        %spectral w-velocity
us = us*(1+1e-9*randn());
u1 = us; u2 = us; unew = us;    %spectral u-velocity for history u-u1-u2
v1 = vs; v2 = vs; vnew = vs;
w1 = ws; w2 = ws; wnew = ws;
u31 = us; u32 = us;             %1/3 space and 2/3 space
v31 = vs; v32 = vs;
w31 = ws; w32 = ws;
% delta1, delta2;                 %correct coefficient
% pc = zeros(nx, ny, nz);         %spectral pressure correct value
% uc = zeros(nx, ny, nz);         %spectral u-velocity correct value
% vc = zeros(nx, ny, nz);         %spectral v-velocity correct value
% wc = zeros(nx, ny, nz);         %spectral w-velocity correct value
[fx, fy, fz] = lamb(us, vs, ws, alpha, beta);         %spectral non-linear term in x-direction
fx(nx/2+1,1,nz/2+1) = fx(nx/2+1,1,nz/2+1)+dpdx*nx*nz;
fx1 = fx; fx2 = fx;             %spectral non-linear term for history fx-fx1-fx2
fy1 = fy; fy2 = fy;
fz1 = fz; fz2 = fz;
div_history = zeros(1,step_num + 2);
dpdx_history = zeros(1,step_num+2);
dpdx_history(1) = dpdx;
%% Main Program
[pc, uc, vc, wc, M] = correct(nx, ny, nz, alpha, beta, re, dt);%compute the correct field
%Main time iteration
for t = start_time:dt:end_time
%     if t > dt*150
%         pause;
%     end
    
    %Mass flux correct
     %t
    if mass_flux ~= 0
        if t == start_time
            temp = dpdx;
            dpdx = dpdx_correct(us, dpdx, zeros(nx, ny, nz), 0, mass_flux, beta);
            dpdx_hist = temp;
        else
            temp = dpdx;
            dpdx = dpdx_correct(us, dpdx, u1, dpdx_hist, mass_flux, beta);
            dpdx_hist = temp;
        end
    end
    [fx,  fy,  fz] = lamb(us, vs, ws, alpha, beta);
    fx(nx/2+1,1,nz/2+1) = fx(nx/2+1,1,nz/2+1)+dpdx*nx*nz;
    
    %Non-linear step
    u31 = 3*us - 3/2*u1 + 1/3*u2 + dt*(3*fx - 3*fx1 + fx2);
    v31 = 3*vs - 3/2*v1 + 1/3*v2 + dt*(3*fy - 3*fy1 + fy2);
    w31 = 3*ws - 3/2*w1 + 1/3*w2 + dt*(3*fz - 3*fz1 + fz2);
    
    %Pressure step
    boundary = 3*Gterm(fy, us, vs, ws, alpha, beta, re) ...
        -3*Gterm(fy1, u1, v1, w1, alpha, beta, re)...
        +  Gterm(fy2, u2, v2, w2, alpha, beta, re);
    boundary = y_trans(boundary, -1);    %transform to (m, y, n) space
    upper_bound = boundary(:,1,:);
    lower_bound = boundary(:,ny,:);
    parfor iter_x = 1:nx
        for iter_z = 1:nz
            ps(iter_x,:,iter_z) = ode2_solve(alpha^2*(iter_x-1-nx/2)^2 + beta^2*(iter_z-1-nz/2)^2,...
                (complex(0,1)*alpha*(iter_x-1-nx/2)*u31(iter_x,:,iter_z)...
                +partial(v31(iter_x,:,iter_z),'y',0)...
                +complex(0,1)*beta*(iter_z-1-nz/2)*w31(iter_x,:,iter_z))/dt,...
                2, [lower_bound(iter_x,1,iter_z) , upper_bound(iter_x,1,iter_z)],...
                's','s');
            u32(iter_x,:,iter_z) = u31(iter_x,:,iter_z) - complex(0,1)*alpha*(iter_x-1-nx/2)*ps(iter_x,:,iter_z)*dt;
            v32(iter_x,:,iter_z) = v31(iter_x,:,iter_z) - partial(ps(iter_x,:,iter_z),'y',0)*dt;
            w32(iter_x,:,iter_z) = w31(iter_x,:,iter_z) - complex(0,1)*beta*(iter_z-1-nz/2)*ps(iter_x,:,iter_z)*dt;
        end
    end
    
    %Viscosity step
    parfor iter_x = 1:nx
        for iter_z = 1:nz
            k = alpha^2*(iter_x-1-nx/2)^2 + beta^2*(iter_z-1-nz/2)^2 + 11/6*re/dt;
            unew(iter_x,:,iter_z) = ode2_solve(k, -re/dt*u32(iter_x,:,iter_z),...
                1, [0,0], 's', 's');
            vnew(iter_x,:,iter_z) = ode2_solve(k, -re/dt*v32(iter_x,:,iter_z),...
                1, [0,0], 's', 's');
            wnew(iter_x,:,iter_z) = ode2_solve(k, -re/dt*w32(iter_x,:,iter_z),...
                1, [0,0], 's', 's');
        end
    end
    
    %Correct
    div_V = globe_partial(unew, 'x', alpha) + globe_partial(vnew, 'y', 0) ...
        + globe_partial(wnew, 'z', beta);
    div_V_p = global_trans(div_V, -1);
    div_V_p = reshape(div_V_p(nx/2+1,:,nz/2+1), 1, ny);
    div_V_p = sum(abs(div_V_p));
    div_history(round((t-start_time)/dt)+1) = div_V_p;
    div_V = y_trans(div_V, -1);
%     if abs(div_V_p)/ny >= 1e-13
%         for iter_x = 1:nx
%             for iter_z = 1:nz
%                 M_temp = reshape(M(iter_x, iter_z, :, :), 2, 2);
%                 div_temp = -1*[div_V(iter_x, 1, iter_z); div_V(iter_x, ny, iter_z)];
%                 if cond(M_temp) > 1e4
% %                     disp('warning');
%                     delta = [div_temp(1)/M_temp(1,1);0];
%                 else
%                     delta = M_temp\div_temp;
%                 end
%                 ps(iter_x,:,iter_z) = ps(iter_x,:,iter_z) + delta(1)*pc(iter_x,:,iter_z)...
%                     + delta(2)*pc(iter_x,:,iter_z).*(-1).^(0:(ny-1));
%                 unew(iter_x,:,iter_z) = unew(iter_x,:,iter_z) + delta(1)*uc(iter_x,:,iter_z)...
%                     + delta(2)*uc(iter_x,:,iter_z).*(-1).^(0:(ny-1));
%                 vnew(iter_x,:,iter_z) = vnew(iter_x,:,iter_z) + delta(1)*vc(iter_x,:,iter_z)...
%                     + delta(2)*vc(iter_x,:,iter_z).*(-1).^(0:(ny-1));
%                 wnew(iter_x,:,iter_z) = wnew(iter_x,:,iter_z) + delta(1)*wc(iter_x,:,iter_z)...
%                     + delta(2)*wc(iter_x,:,iter_z).*(-1).^(0:(ny-1));
%             end
%         end
%     end
    
    p = global_trans(ps, -1);
    u = global_trans(unew, -1);
    v = global_trans(vnew, -1);
    w = global_trans(wnew, -1);
    %save_data(p, u, v, w, t);
    
    %Update
    u2 = u1; u1 = us; us = unew;
    v2 = v1; v1 = vs; vs = vnew;
    w2 = w1; w1 = ws; ws = wnew;
    fx2 = fx1; fx1 = fx;
    fy2 = fy1; fy1 = fy;
    fz2 = fz1; fz1 = fz;
    
    err = u(nx/2, :, nz/2) - u0(nx/2, :, nz/2);
    err = reshape(err, 1, ny);
    dpdx_history(round((t-start_time)/dt)+2) = dpdx;
    subplot(3,1,1)
    plot(y, err);
    subplot(3,1,2)
    plot(dpdx_history);
    subplot(3,1,3)
    plot(div_history);
end