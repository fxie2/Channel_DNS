function y = ysolve(k, r, x, bd1, bd2)
n = length(x) - 1;
dx = zeros(1, n);
for i = 1 :n
    dx(i) = x(i+1) - x(i);
end
h = zeros(1, n+1);
h(1) = dx(1);
for i = 2 : n
    h(i) = (dx(i-1) + dx(i)) / 2;
end
h(n+1) = dx(n);
pfac = zeros(1, n);
cfac = zeros(1, n);
mfac = zeros(1, n);
for i = 1 : n
    pfac(i) = 1 / h(i+1) / (h(i) + h(i+1)) * 2;
    cfac(i) = -2 * (1/h(i+1) + 1/h(i)) / (h(i) + h(i+1)) + k;
    mfac(i) = 2 / h(i) / (h(i) + h(i+1));
end
r1 = r(1);
r(1) = r(1) + bd1 * 2 / (h(1) + h(2));
r(n) = r(n) - bd2 * 2 / (h(n) + h(n+1));
pfac(1) = 2 / h(2) / (h(1) + h(2));
cfac(1) =-2 / h(2) / (h(1) + h(2)) + k;
mfac(1) = 0;
pfac(n) = 0;
cfac(n) = -2 / h(n) / (h(n) + h(n+1)) + k;
mfac(n) = 2 / h(n) / (h(n) + h(n+1));
if (abs(k) < eps())
    pfac(1) = 0;
    cfac(1) = 1;
    mfac(1) = 0;
    r(1) = 0;
    disp('match');
end
M = diag(cfac) + diag(pfac(1:n-1), 1) + diag(mfac(2:n), -1);
if size(r, 1) == 1
    y = M\r.';
%     disp(sum(M*y - r.'));
else
    y = M\r;
%     disp(sum(M*y - r));
end
% if(abs(k) < eps())
%     mean = sum(y.*reshape(dx, size(y)));
%     y = y + (r1 - mean) / (x(n+1) - x(1));
% end
end