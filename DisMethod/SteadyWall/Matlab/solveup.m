global u, v, w, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh;
function [u, v, w, p, pgx, pgz] = solveup(u, v, w, p, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh)
[du, dv, dw] = getvel(u, v, w, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh);
dp = getpre(du, dv, dw, p, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh);
[u, v, w, p, pgx, pgz] = updateup(du, dv, dw, dp, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh);
end

function [du, dv, dw] = getvel(u, v, w, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh)
du = getu(u, v, w, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh);
dv = getv(du, v, w, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh);
dw = getw(du, dv, w, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh);
[du, dv, dw] = finish_vel(du, dv, dw, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh);

du = du + u;
dv = dv + v;
dw = dw + w;
end

function du = getu(u, v, w, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh)
r1 = form_r1(u, v, w, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh);
du = solve_du(u, v, w, r1, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh);
end

function dv = getv(du, v, w, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh)
r2 = form_r2(du, v, w, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh);
dv = solve_dv(du, v, w, r2, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh);
end

function dw = getw(du, dv, w, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh)
r3 = form_r3(du, dv, w, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh);
dw = solve_dw(du, dv, w, r3, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh);
end

function [du, dv, dw] = finish_vel(du, dv, dw, pgx, pgz, dx, dy, dz, h, dy2dy, dy2h, dydy, dyh)

%finish dv
for k = 1:n3
    km = mod(k-2, n3) + 1;
    kp = mod(k, n3) + 1;
    for j = 2:n2
        for i = 1:n1
            
            v1 = (v(i, j, km) + v(i, j, k)) / 2;
            v2 = (v(i, j, kp) + v(i, j, k)) / 2;
            w1 = (dw(i, j+1, k) * dy(