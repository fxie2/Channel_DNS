function [xflow, zflow] = check_flow_rate(u, w, dx, dy, dz)
[n1, n2, n3] = size(u);
n2 = n2 - 2;
xmf = zeros(n1);
zmf = zeros(n3);
for i = 1:n1
    for k = 1:n3
        for j = 1:n2
            xmf(i) = xmf(i) + u(i, j+1, k) * dz * dy(j+1);
        end
    end
end

for k = 1:n3
    for j = 1:n2
        for i = 1:n1
            zmf(k) = zmf(k) + w(i, j+1, k) * dx * dy(j+1);
        end
    end
end

xflow = sum(xmf) / n1;
zflow = sum(zmf) / n3;
end