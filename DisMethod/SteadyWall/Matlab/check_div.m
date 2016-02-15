function divmax = check_div(u, v, w, dx, dy, dz)
divmax = max(max(max(abs(getdiv(u, v, w, dx, dy, dz)))));
end