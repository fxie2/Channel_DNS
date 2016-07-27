%% Channel DNS Subfunction - ode2_solve
%% Purpose
%   Solve the 2nd order diffential equation in the form:
%       f"(y) - k*f(y) = s(y) , y in [-1, 1]
%   with first kind boundary condition:
%       f(-1) = A, f(1) = B
%   or with second kind boundary condition:
%       f'(-1) = A, f'(1) = B
%% Method
%   Use Chebyshev transform to solve the equation.
%   Let
%
%           f(y) = sum fp*Tp(y) , p=0:P
%           s(y) = sum sp*Tp(y) , p=0:P
%
%   The main equation can be written as
%
%   k*c_p-2 * f_p-2          k * e_p+2            k * e_p+4 * f_p+2
%   ---------------- - (1 + ------------ )*f_p + -------------------
%     4*p*(p-1)              2*(p^2-1)                4*p*(p+1)
%
%         c_p-2 * s_p-2      e_p+2 * s_p     e_p+4 * s_p+2
%   = -( ---------------- - ------------- + --------------- )
%            4*p*(p-1)        2*(p^2-1)        4*p*(p+1)
%                                           (by Gottlib & Orszag, 1977)
%   ( e_p = 1 for p<=P, e_p = 0 for p>P)
%   and for boundary conditions
%   1st kind:
%           sum (-1)^p * fp = A , p=0:P
%           sum fp = B , p=0:P
%   2nd kind:
%           sum (-1)^p * fp' = A , p=0:P
%           sum fp' = B , p=0:P
%           cp*fp' = 2 * sum k*fk , k=p+1:P
%% Parameters
%   Input parameters:
%   k ------------------------- the k factor
%   s ------------------------- the right term
%   boundary_type ------------- 1 for 1st, 2 for 2nd.
%   boundary_conditions ------- [A, B] 2 boundary values for (-1, +1)
%   input_type ---------------- 's' for spectral, 'p' for physical
%   output_type --------------- 's' for spectral, 'p' for physical
%   Output parameter:
%   f ------------------------- the solution
%% Author
%   Written by Luohao Wang on 2015-9-10
%   Contact : lh-wang13@mails.tsinghua.edu.cn

%% Code
function f = ode2_solve(k, s, boundary_type, boundary_conditons, input_type, output_type)
P = length(s) - 1;
if input_type == 'p'
    sp = FCT(s, 1);
else
    sp = s;
end
%form the index array
i = zeros(5*P - 3, 1);
p = zeros(5*P - 3, 1);
v = zeros(5*P - 3, 1);
n = 1;
for iter = 1:P-5
    i(n:n+2) = iter*ones(1,3);
    p(n:n+2) = [0,2,4] + iter*ones(1,3);
    v(n) = k / (4*(iter + 1)*iter);
    v(n+1) = -(1 + k / (2*((iter+1)^2 - 1)));
    v(n+2) = k / (4*(iter+1)*(iter+2));
    n = n + 3;
end
for iter = P-4:P-3
    i(n:n+2) = iter*ones(1,3);
    p(n:n+2) = [0,2,4] + iter*ones(1,3);
    v(n) = k / (4*(iter + 1)*iter);
    v(n+1) = -(1 + k / (2*((iter+1)^2 - 1)));
    v(n+2) = 0;
    n = n + 3;
end
v(1) = k / 4;
i(n:n+3) = [P-2,P-2,P-1,P-1];
p(n:n+3) = [P-2,P,P-1,P+1];
v(n) = k/(4*(P-1)*(P-2));
v(n+1) = -1;
v(n+2) = k/(4*P*(P-1));
v(n+3) = -1;
i(3*P-4:4*P-4) = P*ones(1,P+1);
i(4*P-3:5*P-3) = (P+1)*ones(1,P+1);
p(3*P-4:4*P-4) = 1:P+1;
p(4*P-3:5*P-3) = 1:P+1;
if boundary_type == 1
    v(3*P-4:4*P-4) = (-1).^(0:P);
    v(4*P-3:5*P-3) = ones(1,P+1);
elseif boundary_type == 2
    v(3*P-4:4*P-4) = -(0:P).^2.*(-1).^(0:P);
    v(4*P-3:5*P-3) = (0:P).^2;
else
    disp('boundary type error!');
end
M = sparse(i,p,v);
b = zeros(P+1,1);
b(1) = -(sp(1)/4 - sp(3)/6 + sp(5)/24);
for iter = 2:P-5
    pp = iter + 1;
    b(iter) = -(sp(pp-2+1)/(4*pp*(pp-1)) - sp(pp+1)/(2*(pp^2 - 1)) + sp(pp+3)/(4*pp*(pp+1)));
end
b(P-4) = -(sp(P-4)/(4*(P-3)*(P-4)) - sp(P-2)/(2*((P-3)^2 - 1)));
b(P-3) = -(sp(P-3)/(4*(P-2)*(P-3)) - sp(P-1)/(2*((P-2)^2 - 1)));
b(P-2) = -sp(P-2)/(4*(P-1)*(P-2));
b(P-1) = -sp(P-1)/(4*(P-1)*P);
b(P) = boundary_conditons(1);
b(P+1) = boundary_conditons(2);
if abs(k) <= eps() && boundary_type == 2
    temp_ans = zeros(P+1,1);
    temp_ans(2:P+1) = M(2:P+1,2:P+1)\b(2:P+1);
else
    temp_ans = M\b;
end
if output_type == 'p'
    f = FCT(temp_ans,-1);
else
    f = temp_ans;
end
end