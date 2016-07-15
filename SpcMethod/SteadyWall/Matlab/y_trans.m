function f = y_trans(u, direction)
[nx, ~, nz] = size(u);
f = zeros(size(u));
for iter_x = 1:nx
    for iter_z = 1:nz
        f(iter_x,:,iter_z) = FCT(u(iter_x,:,iter_z), direction);
    end
end
end