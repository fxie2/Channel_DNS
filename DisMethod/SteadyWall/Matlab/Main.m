%% Channel DNS

%% Init parameters
n1 = 16;
n2 = 128;
n3 = 32;
turb_intensity = 0.01;
re = 1000;
dt = 0.01;
start_time = 0;
end_time = 100;


%% Mesh
lx = 2 * pi;
lz = 0.5 * pi;
ly = 2;
y = linspace(-1, 1, n2 + 1);
x = linspace(0, lx, n1 + 1);
z = linspace(0, lz, n3 + 1);
dx = lx / n1;
dz = lz / n3;
dy = zeros(1, n2 + 2);
dy(2:n2+1) = y(2:n2+1) - y(1:n2);
h = zeros(1, n2 + 1);
h(1:n2 + 1) = (dy(1:n2+1) + dy(2:n2+2)) / 2;
dy2dy = zeros(3, n2);
dy2dy(1, :) = 1 ./ h(1:n2) ./ dy(2:n2+1);
dy2dy(2, :) = -(1 ./ dy(2:n2+1) + 1 ./ dy(1:n2)) ./ h(1:n2);
dy2dy(3, :) = 1 ./ h(1:n2) ./ dy(1:n2);
dydy = zeros(3, n2);
dydy(1, :) = 1 ./ dy(2:n2+1) ./ h(1:n2);
dydy(2, :) = -(1 ./ dy(2:n2+1) + 1 / dy(1:n2)) ./ h(1:n2);
dydy(3, :) = 1 ./ dy(1:n2) ./ h(1:n2);
dy2h = zeros(3, n2);
dy2h(1, :) = 1 ./ h(2:n2+1) / dy(2:n2+1);
dy2h(2, :) = -(h(2:n2+1) + h(1:n2)) ./ h(2:n2+1) ./ h(1:n2) ./ dy(2:n2+1);
dy2h(3, :) = 1 ./ h(1:n2) ./ dy(2:n2+1);
dyh = zeros(3, n2);
dyh(1, :) = h(1:n2) .* h(2:n2+1) ./ (h(1:n2) + h(2:n2+1)) ./ h(2:n2+1) ./ h(2:n2+1);
dyh(2, :) = h(1:n2) .* h(2:n2+1) ./ (h(1:n2) + h(2:n2+1)) .* (1 ./ h(1:n2) ./ h(1:n2) - 1 ./ h(2:n2+1) .^ 2);
dyh(3, :) = h(1:n2) .* h(2:n2+1) ./ (h(1:n2) + h(2:n2+1)) ./ h(1:n2) .^ 2 .* -1.0;

%% Field
p = zeros(n1, n2, n3);
pgx = -2 / re;
pgz = 0;
div = zeros(n1, n2, n3);

%% Init Field
u = rand(n1, n2 + 2, n3);
v = rand(n1, n2 + 1, n3);
w = rand(n1, n2 + 2, n3);
u = u - 0.5;
v = v - 0.5;
w = w - 0.5;
u(:, 1, :) = 0;
v(:, 1, :) = 0;
w(:, 1, :) = 0;
u(:, n2 + 2, :) = 0;
v(:, n2 + 1, :) = 0;
w(:, n2 + 2, :) = 0;
for i = 1:n1
    v1m = 0;
    s1 = 0;
    for k = 1:n3
       for j = 1:n2
           s1 = s1 + dy(j+1) .* dz;
           v1m = v1m + u(i, j+1, k) .* dy(j+1) .* dz;
       end
    end
    v1m = v1m / s1;
    u(i, 2:n2+1, :) = u(i, 2:n2+1, :) - v1m;
end

for j = 2:n2
    v2m = 0;
    s2 = 0;
    for k = 1:n3
        for i = 1:n1
            s2 = s2 + dx * dz;
            v2m = v2m + v(i, j, k) * dx * dz;
        end
    end
    v2m = v2m / s2;
    v(:, j, :) = v(:, j, :) - v2m;
end

for k = 1:n3
    v3m = 0;
    s3 = 0;
    for j = 1:n2
        for i = 1:n1
            s3 = s3 + dy(j+1) * dx;
            v3m = v3m + w(i, j+1, k) * dy(j+1) * dx;
        end
    end
    v3m = v3m / s3;
    w(:, 2:n2+1, k) = w(:, 2:n2+1, k) - v3m;
end

u = u * turb_intensity * 2;
v = v * turb_intensity * 2;
w = w * turb_intensity * 2;

for j = 2:n2+1
    yh = (y(j) + y(j - 1)) / 2;
    u(:, j, :) = u(:, j, :) + 1 - yh .^ 2;
end

rflow = 4.0 / 3.0 * lz;
[xflow, zflow] = check_flow_rate(u, w, dx, dy, dz);
u = rflow / xflow * u;
w(:, 2:n2+1, :) = w(:, 2:n2+1, :) - zflow / lx / ly;

t = 0;
while(t < end_time)
    [u, v, w, p, pgx, pgz] = solveup(u, v, w, p, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh);
    disp(check_div(u, v, w, dx, dy, dz));
end