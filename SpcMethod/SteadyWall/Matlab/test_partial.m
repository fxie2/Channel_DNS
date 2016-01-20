y = cos((0:64)/64*pi);
a = exp(y);
as = FCT(a, 1);
das = partial(as, 'y', 0);
da = FCT(das, -1);
err = da - a;
plot(y, real(err));