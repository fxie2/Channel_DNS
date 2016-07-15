function div = getdiv(u, v, w, dx, dy, dz)
[n1, n2, n3] = size(u);
n2 = n2 - 2;

div = zeros(n1, n2, n3);

for k = 1:n3
    kp = mod(k, n3) + 1;
    for j = 1:n2
        jp = j+1;
        for i = 1:n1
            ip = mod(i, n1) + 1;
            div(i, j, k) = (u(ip, j+1, k) - u(i, j+1, k)) / dx ...
                         + (v(i, jp, k) - v(i, j, k)) / dy(j) ...
                         + (w(i, j+1, kp) - w(i, j+1, k)) / dz;
        end
    end
end
end