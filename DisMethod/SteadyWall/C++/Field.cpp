#include "Field.h"
#include "Parameter.h"
#include "Mesh.h"
#include "Math.h"
#include <memory>

Field u, v, w, p;
Field du, dv, dw, dp;
Field divs_re, divs_im, dps_re, dps_im;
Field r;
Field &r1 = r, &r2 = r, &r3 = r, &rp = r;
Field dive;

double pgx, pgz;
double xflow, zflow;
double divmax;
double cflmax;

// Supplimental vector
// Factor vector for td & ctd solver
double *nmfac, *ncfac, *npfac;
double *pmfac_re, *pcfac_re, *ppfac_re;
double *pmfac_im, *pcfac_im, *ppfac_im;

void alloc_field()
{
	u.reAlloc(n1, n2 + 1, n3, 1, 0, 1);
	v.reAlloc(n1, n2 + 1, n3);
	w.reAlloc(n1, n2 + 1, n3, 1, 0, 1);
	p.reAlloc(n1, n2, n3);
	dive.reAlloc(n1, n2, n3);
	du.reAlloc(n1, n2 + 1, n3, 1, 0, 1);
	dv.reAlloc(n1, n2 + 1, n3);
	dw.reAlloc(n1, n2 + 1, n3, 1, 0, 1);
	dp.reAlloc(n1, n2, n3);
	dps_re.reAlloc(n1, n2, n3);
	dps_im.reAlloc(n1, n2, n3);
	divs_re.reAlloc(n1, n2, n3);
	divs_im.reAlloc(n1, n2, n3);
	r.reAlloc(n1, n2, n3);
	int sizeMax = (n1 > n2) ? n1 : n2;
	sizeMax = (sizeMax > n3) ? sizeMax : n3;
	nmfac = new double[sizeMax];
	ncfac = new double[sizeMax];
	npfac = new double[sizeMax];
	pmfac_re = new double[sizeMax];
	pmfac_im = new double[sizeMax];
	pcfac_re = new double[sizeMax];
	pcfac_im = new double[sizeMax];
	ppfac_re = new double[sizeMax];
	ppfac_im = new double[sizeMax];
}

void dealloc_field()
{
	delete[] nmfac, ncfac, npfac;
	delete[] pmfac_im, pmfac_re, pcfac_im, pcfac_re, ppfac_im, ppfac_re;
}

void iniup()
{
	double v1m, v2m, v3m;	//mean flow
	double s1, s2, s3;		//slice area
	double yh;
	double rflow;

	//uniform distribute in [-0.5, 0.5)
	u.setRand();
	v.setRand();
	w.setRand();

	u -= 0.5;
	v -= 0.5;
	w -= 0.5;

	//impose zero velocity at boundary
	for (size_t i = 1; i < n1 + 1; i++)
	{
		for (size_t k = 1; k < n3 + 1; k++)
		{
			u(i, 0, k) = 0;
			v(i, 1, k) = 0;
			w(i, 0, k) = 0;
			u(i, n2 + 1, k) = 0;
			v(i, n2 + 1, k) = 0;
			w(i, n2 + 1, k) = 0;
		}
	}

	//eliminate mean quantities of random fluctuations
	//u direction
	for (size_t i = 1; i < n1+1; i++)
	{
		v1m = 0;
		s1 = 0;
		for (size_t j = 1; j < n2+1; j++)
		{
			for (size_t k = 1; k < n3+1; k++)
			{
				s1 += dy[j] * dz;
				v1m += u(i, j, k) * dy[j] * dz;
			}
		}
		v1m /= s1;
		for (size_t j = 1; j < n2+1; j++)
		{
			for (size_t k = 1; k < n3+1; k++)
			{
				u(i, j, k) -= v1m;
			}
		}
	}

	//v direction
	for (size_t j = 2; j < n2+1; j++)
	{
		v2m = 0;
		s2 = 0;
		for (size_t i = 1; i < n1+1; i++)
		{
			for (size_t k = 1; k < n3+1; k++)
			{
				s2 += dx * dz;
				v2m += v(i, j, k) * dx * dz;
			}
		}
		v2m /= s2;
		for (size_t i = 1; i < n1+1; i++)
		{
			for (size_t k = 1; k < n3+1; k++)
			{
				v(i, j, k) -= v2m;
			}
		}
	}

	//w direction
	for (size_t k = 1; k < n3+1; k++)
	{
		v3m = 0;
		s3 = 0;
		for (size_t i = 1; i < n1+1; i++)
		{
			for (size_t j = 1; j < n2+1; j++)
			{
				s3 += dy[j] * dx;
				v3m += w(i, j, k) * dy[j] * dx;
			}
		}
		v3m /= s3;
		for (size_t i = 1; i < n1+1; i++)
		{
			for (size_t j = 1; j < n2+1; j++)
			{
				w(i, j, k) -= v3m;
			}
		}
	}

	u *= init_turb_intensity * 2;
	v *= init_turb_intensity * 2;
	w *= init_turb_intensity * 2;

	//impose laminar velocity profiels in u velocities
	for (size_t j = 1; j < n2+1; j++)
	{
		yh = (y[j] + y[jm[j]]) / 2;
		for (size_t i = 1; i < n1+1; i++)
		{
			for (size_t k = 1; k < n3+1; k++)
			{
				u(i, j, k) += (1 - yh * yh);
			}
		}
	}

	rflow = 4.0 / 3.0 * lz; //integral of u along y
	check_flow_rate();
	u *= rflow / xflow;
	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t j = 1; j < n2 + 1; j++)
			for (size_t k = 1; k < n3 + 1; k++)
				w(i, j, k) -= zflow / lx / ly;

	p = 0;

	pgx = -2 / re;
	pgz = 0;
}

void check_flow_rate()
{
	double *xmf = new double[n1 + 1];
	double *zmf = new double[n3 + 1];

	for (size_t i = 1; i < n1+1; i++)
	{
		xmf[i] = 0;
		for (size_t j = 1; j < n2 + 1; j++)
			for (size_t k = 1; k < n3 + 1; k++)
				xmf[i] += u(i, j, k) * dz * dy[j];
	}

	for (size_t k = 0; k < n3+1; k++)
	{
		zmf[k] = 0;
		for (size_t i = 1; i < n1 + 1; i++)
			for (size_t j = 1; j < n2 + 1; j++)
				zmf[k] += w(i, j, k) * dx * dy[j];
	}

	xflow = 0;
	for (size_t i = 1; i < n1 + 1; i++)
		xflow += xmf[i];
	xflow /= n1;

	zflow = 0;
	for (size_t k = 1; k < n3 + 1; k++)
		zflow += zmf[k];
	zflow /= n3;

	delete[] xmf, zmf;
}

void check_div()
{
	getdiv(u, v, w, dive, t);
	divmax = 0;
	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t j = 1; j < n2 + 1; j++)
			for (size_t k = 1; k < n3 + 1; k++)
				divmax = (divmax > abs(dive(i, j, k))) ? divmax : abs(dive(i, j, k));
}

void getdiv(Field& uu, Field& vv, Field& ww, Field& ddiv, double tt)
{
	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t j = 1; j < n2 + 1; j++)
			for (size_t k = 1; k < n3 + 1; k++)
				ddiv(i, j, k) = (uu(ip[i], j, k) - uu(i, j, k)) / dx
				+ (vv(i, jp[j], k) - vv(i, j, k)) / dy[j]
				+ (ww(i, j, kp[k]) - ww(i, j, k)) / dz;
}

void solveup()
{
	getvel();
	getpre();
	update_up();
}

void getvel()
{
	update_bc();

	getu();
	getv();
	getw();
	finish_vel();

	du += u;
	dv += v;
	dw += w;
}

void update_bc()
{
	du = 0;
	dv = 0;
	dw = 0;
}

void getpre()
{
	form_rp();
	solve_dp();
}

void getu()
{
	form_r1();
	solve_du();
}

void getv()
{
	form_r2();
	solve_dv();
}

void getw()
{
	form_r3();
	solve_dw();
}

void finish_vel()
{
	double v1, v2, w1, w2, dwvdz, dwvdy;
	double u1, u2, w_up, w_mi, w_dn, duvdy, duwdz, duwdy;

	//finish dv
	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t j = 2; j < n2 + 1; j++)
			for (size_t k = 1; k < n3+1; k++)
			{
				//d_w^np_v^n_dz
				v1 = (v(i, j, km[k]) + v(i, j, k)) / 2;
				v2 = (v(i, j, kp[k]) + v(i, j, k)) / 2;
				w1 = (dw(i, j, k) * dy[jm[j]] + dw(i, jm[j], k) * dy[j]) / 2 / h[j];
				w2 = (dw(i, j, kp[k]) * dy[jm[j]] + dw(i, jm[j], kp[k]) * dy[j]) / 2 / h[j];
				dwvdz = (w2 * v2 - w1 * v1) / dz;

				dv(i, j, k) = dv(i, j, k) - dt * dwvdz / 2;
			}

	//finish du
	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t j = 1; j < n2 + 1; j++)
			for (size_t k = 1; k < n3 + 1; k++)
			{
				//d_u^n_v^n_dy
				u1 = (u(i, j, k) * dy[jm[j]] + u(i, jm[j], k) * dy[j]) / 2 / h[j];
				u2 = (u(i, j, k) * dy[jp[j]] + u(i, jp[j], k) * dy[j]) / 2 / h[jp[j]];
				v1 = (dv(im[i], j, k) + dv(i, j, k)) / 2;
				v2 = (dv(im[i], jp[j], k) + dv(i, jp[j], k)) / 2;
				duvdy = (u2 * v2 - u1 * v1) / dy[j];

				//d_u^n_w^n_dz
				u1 = (u(i, j, k) + u(i, j, km[k])) / 2;
				u2 = (u(i, j, kp[k]) + u(i, j, k)) / 2;
				w1 = (dw(im[i], j, k) + dw(i, j, k)) / 2;
				w2 = (dw(i, j, kp[k]) + dw(im[i], j, kp[k])) / 2;
				duwdz = (u2 * w2 - u1 * w1) / dz;

				du(i, j, k) = du(i, j, k) - dt * duvdy / 2 - dt * duwdz / 2;
			}
}

void form_r1()
{
	double viscos, cross, nonlin, pressg;
	double dudx[4], dudz[4], dudy;
	double duudx, duvdy, duwdz, duudy, duwdy;
	double u1, u2, v1, v2, w1, w2, w_up, w_mi, w_dn;
	double xc, yc, zc;
	double vis1, vis2, vis3, u11, u22, u33, fac1, fac2, fac3;
	double bc_dn, bc_up, bcond;

	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t j = 1; j < n2 + 1; j++)
			for (size_t k = 1; k < n3 + 1; k++)
			{
				//viscos term
				viscos = (u(ip[i], j, k) - 2 * u(i, j, k) + u(im[i], j, k)) / dx / dx
					+ (dy2h1[j] * u(i, jp[j], k) + dy2h2[j] * u(i, j, k) + dy2h3[j] * u(i, jm[j], k))
					+ (u(i, j, kp[k]) - 2 * u(i, j, k) + u(i, j, km[k])) / dz / dz;
				//form nth-time-step viscos term
				viscos = viscos / re * dt;

				//nolinear term
				//d_u^n_u^n_dx
				u1 = (u(im[i], j, k) + u(i, j, k)) / 2;
				u2 = (u(ip[i], j, k) + u(i, j, k)) / 2;
				duudx = (u2 * u2 - u1 * u1) / dx;

				//d_u^n_v^n_dy
				u1 = (u(i, j, k) * dy[jm[j]] + u(i, jm[j], k) * dy[j]) / 2 / h[j];
				u2 = (u(i, j, k) * dy[jp[j]] + u(i, jp[j], k) * dy[j]) / 2 / h[jp[j]];
				v1 = (v(im[i], j, k) + v(i, jp[j], k)) / 2;
				v2 = (v(im[i], jp[j], k) + v(i, jp[j], k)) / 2;
				duvdy = (u2 * v2 - u1 * v1) / dy[j];

				//d_u^n_w^n_dz
				u1 = (u(i, j, k) + u(i, j, km[k])) / 2;
				u2 = (u(i, j, kp[k]) + u(i, j, k)) / 2;
				w1 = (w(im[i], j, k) + w(i, j, k)) / 2;
				w2 = (w(i, j, kp[k]) + w(im[i], j, kp[k])) / 2;
				duwdz = (u2 * w2 - u1 * w1) / dz;

				nonlin = duudx + duvdy + duwdz;

				//boundary term
				bc_dn = dy2h3[1] * du(i, 0, k) / re / 2
					+ du(i, 0, k) * (v(i, 1, k) + v(im[i], 1, k)) / 2 / dy[1] / 2
					+ u(i, 0, k) * (dv(i, 1, k) + dv(im[i], 1, k)) / 2 / dy[1] / 2;

				bc_up = dy2h1[n2] * du(i, n2 + 1, k) / re / 2
					- du(i, n2 + 1, k) * (v(i, n2 + 1, k) + v(im[i], n2 + 1, k)) / 2 / dy[n2] / 2
					- u(i, n2 + 1, k) * (dv(i, n2 + 1, k) + dv(im[i], n2 + 1, k)) / 2 / dy[n2] / 2;

				bcond = ((n2 - j + 1) - (n2 - j + 1) % n2) / (double)n2 * bc_dn + (j - j%n2) / (double)n2 * bc_up;

				//pressure term
				pressg = (p(i, j, k) - p(im[i], j, k)) / dx + pgx;

				//form r1 term
				r1(i, j, k) = viscos - nonlin * dt - pressg * dt - bcond * dt;
			}
}

void form_r2()
{
	double viscos, cross, nonlin, pressg;
	double dvdx[4], dvdz[4], dvdy;
	double duvdx, dvvdy, dwvdz, duvdy, dwvdy;
	double u1, u2, v1, v2, w1, w2;
	double m21u;
	double bc_dn, bc_up, bcond;

	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t j = 2; j < n2 + 1; j++)
			for (size_t k = 1; k < n3 + 1; k++)
			{
				//viscos term
				viscos = (v(ip[i], j, k) - 2 * v(i, j, k) + v(im[i], j, k)) / dx / dx
					+ (dy2dy1[j] * v(i, jp[j], k) + dy2dy2[j] * v(i, j, k) + dy2dy3[j] * v(i, jm[j], k))
					+ (v(i, j, kp[k]) - 2 * v(i, j, k) + v(i, j, km[k])) / dz / dz;

				viscos = viscos / re * dt;

				//nonlinear term
				//d_u^n_v^n_dx
				u1 = (u(i, jm[j], k) * dy[j] + u(i, j, k) * dy[jm[j]]) / 2 / h[j];
				u2 = (u(ip[i], jm[j], k) * dy[j] + u(ip[i], j, k) * dy[jm[j]]) / 2 / h[j];
				v1 = (v(i, j, k) + v(im[i], j, k)) / 2;
				v2 = (v(i, j, k) + v(ip[i], j, k)) / 2;
				duvdx = (u2 * v2 - u1 * v1) / dx;

				//d_v^n_v^n_dy & d_u^n_v^n_dy
				u1 = (u(i, jm[j], k) + u(ip[i], jm[j], k)) / 2;
				u2 = (u(i, j, k) + u(ip[i], j, k)) / 2;
				v1 = (v(i, jm[j], k) + v(i, j, k)) / 2;
				v2 = (v(i, jp[j], k) + v(i, j, k)) / 2;
				dvvdy = (v2 * v2 - v1 * v1) / h[j];

				//d_w^n_v^n_dz
				v1 = (v(i, j, km[k]) + v(i, j, k)) / 2;
				v2 = (v(i, j, kp[k]) + v(i, j, k)) / 2;
				w1 = (w(i, j, k) * dy[jm[j]] + w(i, jm[j], k) * dy[j]) / 2 / h[j];
				w2 = (w(i, j, kp[k]) * dy[jm[j]] + w(i, jm[j], kp[k]) * dy[j]) / 2 / h[j];
				dwvdz = (w2 * v2 - w1 * v1) / dz;

				nonlin = duvdx + dvvdy + dwvdz;

				//pressure term
				pressg = (p(i, j, k) - p(i, jm[j], k)) / h[j];

				//boundary term
				bc_dn = dv(i, 1, k) * dy2dy3[2] / re / 2
					+ dv(i, 1, k) * (v(i, 1, k) + v(i, 2, k)) / 2 / h[2] / 2;
				bc_up = dv(i, n2 + 1, k) * dy2dy1[n2] / re / 2
					- dv(i, n2 + 1, k) * (v(i, n2 + 1, k) + v(i, n2, k)) / 2 / h[n2] / 2;

				bcond = ((n2 - j + 1) - (n2 - j + 1) % n2) / (double)n2 * bc_dn
					+ (j - j % n2) / (double)n2 * bc_up;

				//m21 term

				//d_v^n_u^np_dx
				u1 = (dy[j] * du(i, jm[j], k) + dy[jm[j]] * du(i, j, k)) / 2 / h[j];
				u2 = (dy[j] * du(ip[i], jm[j], k) + dy[jm[j]] * du(ip[i], j, k)) / 2 / h[j];
				v1 = (v(i, j, k) + v(im[i], j, k)) / 2;
				v2 = (v(ip[i], j, k) + v(i, j, k)) / 2;

				duvdx = (u2 * v2 - u1 * v1) / dx;

				m21u = duvdx / 2;

				//form r2 term
				r2(i, j, k) = viscos
					-nonlin * dt - pressg * dt - bcond * dt - m21u * dt;
			}
}

void form_r3()
{
	double viscos, cross, nonlin, pressg;
	double dwdx[4], dwdz[4], dwdy;
	double duwdx, dvwdy, dwwdz, duwdy, dwwdy;
	double u1, u2, v1, v2, w1, w2, u_up, u_mi, u_dn;
	double m31u, m32v;
	double bc_dn, bc_up, bcond;

	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t j = 1; j < n2 + 1; j++)
			for (size_t k = 1; k < n3 + 1; k++)
			{
				//viscos term
				viscos = (w(ip[i], j, k) - 2 * w(i, j, k) + w(im[i], j, k)) / dx / dx
					+ (dy2h1[j] * w(i, jp[j], k) + dy2h2[j] * w(i, j, k) + dy2h3[j] * w(i, jm[j], k))
					+ (w(i, j, kp[k]) - 2 * w(i, j, k) + w(i, j, km[k])) / dz / dz;

				//form nth - time - step viscos term
				viscos = viscos / re * dt;

				//nonlinear term
				//d_u^n_w^n_dx
				u1 = (u(i, j, kp[k]) + u(i, j, k)) / 2;
				u2 = (u(ip[i], j, kp[k]) + u(ip[i], j, k)) / 2;
				w1 = (w(im[i], j, k) + w(i, j, k)) / 2;
				w2 = (w(ip[i], j, k) + w(i, j, k)) / 2;
				duwdx = (u2 * w2 - u1 * w1) / dx;

				//d_v^n_w^n_dy & d_w^n_w^n_dy
				v1 = (v(i, j, k) + v(i, j, km[k])) / 2;
				v2 = (v(i, jp[j], k) + v(i, jp[j], km[k])) / 2;
				w1 = (w(i, j, k) * dy[jm[j]] + w(i, jm[j], k) * dy[j]) / 2 / h[j];
				w2 = (w(i, jp[j], k) * dy[j] + w(i, j, k) * dy[jp[j]]) / 2 / h[jp[j]];
				dvwdy = (v2 * w2 - v1 * w1) / dy[j];

				//d_w^n_w^n_dz
				w1 = (w(i, j, km[k]) + w(i, j, k)) / 2;
				w2 = (w(i, j, kp[k]) + w(i, j, k)) / 2;
				dwwdz = (w2 * w2 - w1 * w1) / dz;

				nonlin = duwdx + dvwdy + dwwdz;

				//pressure term
				pressg = (p(i, j, k) - p(i, j, km[k])) / dz + pgz;

				//boundary term
				bc_dn = dw(i, 0, k) * dy2h3[1] / re / 2
					+ dw(i, 0, k) * (v(i, 1, k) + v(i, 1, km[k])) / 2 / dy[1] / 2
					+ w(i, 0, k) * (dv(i, 1, k) + dv(i, 1, km[k])) / 2 / dy[1] / 2;

				bc_up = dw(i, n2 + 1, k) * dy2h1[n2] / re / 2
					- dw(i, n2 + 1, k) * (v(i, n2 + 1, k) + v(i, n2 + 1, km[k])) / 2 / dy[n2] / 2
					- w(i, n2 + 1, k) * (dv(i, n2 + 1, k) + dv(i, n2 + 1, km[k])) / 2 / dy[n2] / 2;

				bcond = ((n2 - j + 1) - (n2 - j + 1) % n2) / double(n2) * bc_dn
					+ (j - j % n2) / double(n2) * bc_up;

				//m31 term

				//d_w^n_u^np_dx
				u1 = (du(i, j, kp[k]) + du(i, j, k)) / 2;
				u2 = (du(ip[i], j, kp[k]) + du(ip[i], j, k)) / 2;
				w1 = (w(im[i], j, k) + w(i, j, k)) / 2;
				w2 = (w(i, j, k) + w(ip[i], j, k)) / 2;

				duwdx = (u2 * w2 - u1 * w1) / dx;

				m31u = duwdx / 2;

				//m32 term

				//d_w^n_v^np_dy
				v1 = (dv(i, j, k) + dv(i, j, km[k])) / 2;
				v2 = (dv(i, jp[j], k) + dv(i, jp[j], km[k])) / 2;
				w1 = (w(i, j, k) * dy[jm[j]] + w(i, jm[j], k) * dy[j]) / 2 / h[j];
				w2 = (w(i, jp[j], k) * dy[j] + w(i, j, k) * dy[jp[j]]) / 2 / h[jp[j]];

				dvwdy = (v2 * w2 - v1 * w1) / h[j];

				m32v = dvwdy / 2;

				//form r3 term
				r3(i, j, k) = viscos
					- nonlin * dt - pressg * dt - bcond * dt - m31u * dt - m32v * dt;
			}

}

void solve_du()
{
	double u1, u2, v1, v2, w1, w2, phi, w_up, w_mi, w_dn;
	double visnm, visnc, visnp, visfac;
	double maxdu;
	double hh[34];

	for (size_t i = 0; i < 34; i++)
	{
		hh[i] = h[i];
	}

	//solve in y direction
	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t k = 1; k < n3 + 1; k++)
		{
			for (size_t j = 1; j < n2 + 1; j++)
			{
				//d_v^n_u^np_dy
				v1 = (v(im[i], j, k) + v(i, j, k)) / 2;
				v2 = (v(im[i], jp[j], k) + v(i, jp[j], k)) / 2;
				npfac[j] = v2 / 2 / h[jp[j]];
				ncfac[j] = (v2 * dy[jp[j]] / 2 / h[jp[j]] - v1 * dy[jm[j]]) / 2 / h[j] / dy[j];
				nmfac[j] = -v1 / 2 / h[j];

				//viscos term
				visnp = dy2h1[j];
				visnc = dy2h2[j];
				visnm = dy2h3[j];

				npfac[j] -= visnp / re;
				ncfac[j] -= visnc / re;
				nmfac[j] -= visnm / re;

				npfac[j] = npfac[j] / 2 * dt;
				ncfac[j] = ncfac[j] / 2 * dt + 1;
				nmfac[j] = nmfac[j] / 2 * dt;
			}
			tdma(nmfac, ncfac, npfac, &du(i, 1, k), &r1(i, 1, k), n3, n3, n2);
		}

	//solve in x direction
	for (size_t j = 1; j < n2 + 1; j++)
		for (size_t k = 1; k < n3 + 1; k++)
		{
			for (size_t i = 1; i < n1 + 1; i++)
			{
				//d_u^n_u^np_dx
				u1 = (u(im[i], j, k) + u(i, j, k)) / 2;
				u2 = (u(ip[i], j, k) + u(i, j, k)) / 2;
				npfac[i] = u2 / 2 / dx * dt;
				ncfac[i] = (u2 - u1) / 2 / dx * dt;
				nmfac[i] = -u1 / 2 / dx * dt;

				//viscos term
				npfac[i] -= 1 / dx / dx / 2 / re * dt;
				ncfac[i] += 1 / dx / dx / re * dt;
				nmfac[i] -= 1 / dx / dx / 2 / re * dt;

				ncfac[i] += 1;
			}
			ctdma(nmfac, ncfac, npfac, &du(1, j, k), &du(1, j, k), (n2 + 2)*(n3), (n2 + 2)*(n3), n1);
		}

	//solve in z direction
	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t j = 1; j < n2 + 1; j++)
		{
			for (size_t k = 1; k < n3 + 1; k++)
			{
				//d_w^n_u^np_dz
				w1 = (w(im[i], j, k) + w(i, j, k)) / 2;
				w2 = (w(i, j, kp[k]) + w(im[i], j, kp[k])) / 2;
				npfac[k] = w2 / 2 / dz;
				ncfac[k] = (w2 - w1) / 2 / dz;
				nmfac[k] = -w1 / 2 / dz;

				//viscos term
				npfac[k] -= 1 / dz / dz / re;
				ncfac[k] += 2 / dz / dz / re;
				nmfac[k] -= 1 / dz / dz / re;

				npfac[k] = npfac[k] / 2 * dt;
				ncfac[k] = ncfac[k] / 2 * dt + 1;
				nmfac[k] = nmfac[k] / 2 * dt;
			}
			ctdma(nmfac, ncfac, npfac, &du(i, j, 1), &du(i, j, 1), 1, 1, n3);
		}
}

void solve_dv()
{
	double u1, u2, v1, v2, w1, w2, phi;
	double visnm, visnc, visnp, visfac;

	//solve in y direction
	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t k = 1; k < n3 + 1; k++)
		{
			for (size_t j = 2; j < n2 + 1; j++)
			{
				//d_v^n_v^np_dy
				v1 = (v(i, jm[j], k) + v(i, j, k)) / 2;
				v2 = (v(i, jp[j], k) + v(i, j, k)) / 2;
				npfac[j] = v2 / 2 / h[j];
				ncfac[j] = (v2 - v1) / 2 / h[j];
				nmfac[j] = -v1 / 2 / h[j];

				//viscos term
				visfac = 1. / 2 / re;
				visnp = visfac * dy2dy1[j];
				visnc = visfac * dy2dy2[j];
				visnm = visfac * dy2dy3[j];

				npfac[j] -= visnp;
				ncfac[j] -= visnc;
				nmfac[j] -= visnm;

				npfac[j] = npfac[j] * dt;
				ncfac[j] = ncfac[j] * dt + 1;
				nmfac[j] = nmfac[j] * dt;
			}
			tdma(&nmfac[1], &ncfac[1], &npfac[1], &dv(i, 2, k), &r2(i, 2, k), n3, n3, n2 - 1);
		}

	//solve in x direction
	for (size_t j = 2; j < n2 + 1; j++)
		for (size_t k = 1; k < n3 + 1; k++)
		{
			for (size_t i = 1; i < n1 + 1; i++)
			{
				//d_u^n_v^np_dx
				u1 = (u(i, jm[j], k) * dy[j] + u(i, j, k) * dy[jm[j]]) / 2 / h[j];
				u2 = (u(ip[i], jm[j], k) * dy[j] + u(ip[i], j, k) * dy[jm[j]]) / 2 / h[j];
				npfac[i] = u2 / 2 / dx;
				ncfac[i] = (u2 - u1) / 2 / dx;
				nmfac[i] = -u1 / 2 / dx;

				//viscos term
				npfac[i] -= 1 / dx / dx / re;
				ncfac[i] += 2 / dx / dx / re;
				nmfac[i] -= 1 / dx / dx / re;

				npfac[i] = npfac[i] / 2 * dt;
				ncfac[i] = ncfac[i] / 2 * dt + 1;
				nmfac[i] = nmfac[i] / 2 * dt;
			}
			ctdma(nmfac, ncfac, npfac, &dv(1, j, k), &dv(1, j, k), (n2 + 1)*n3, (n2 + 1)*n3, n1);
		}

	//solve in z direction
	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t j = 2; j < n2 + 1; j++)
		{
			for (size_t k = 1; k < n3 + 1; k++)
			{
				//d_w^n_v^np_dz
				w1 = (w(i, j, k) * dy[jm[j]] + w(i, jm[j], k) * dy[j]) / 2 / h[j];
				w2 = (w(i, j, kp[k]) * dy[jm[j]] + w(i, jm[j], kp[k]) * dy[j]) / 2 / h[j];
				npfac[k] = w2 / 2 / dz;
				ncfac[k] = (w2 - w1) / 2 / dz;
				nmfac[k] = -w1 / 2 / dz;

				//viscos term
				npfac[k] -= 1 / dz / dz / re;
				ncfac[k] += 2 / dz / dz / re;
				nmfac[k] -= 1 / dz / dz / re;

				npfac[k] = npfac[k] / 2 * dt;
				ncfac[k] = ncfac[k] / 2 * dt + 1;
				nmfac[k] = nmfac[k] / 2 * dt;
			}
			ctdma(nmfac, ncfac, npfac, &dv(i, j, 1), &dv(i, j, 1), 1, 1, n3);
		}
}

void solve_dw()
{
	double u1, u2, v1, v2, w1, w2, phi, u_up, u_mi, u_dn;
	double visnm, visnc, visnp, visfac;

	//solve in y direction
	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t k = 1; k < n3 + 1; k++)
		{
			for (size_t j = 1; j < n2 + 1; j++)
			{
				//d_v^n_w^np_dy
				v1 = (v(i, j, k) + v(i, j, km[k])) / 2;
				v2 = (v(i, jp[j], k) + v(i, jp[j], km[k])) / 2;
				npfac[j] = v2 / 2 / h[jp[j]];
				ncfac[j] = (v2 * dy[jp[j]] / h[jp[j]] - v1 * dy[jm[j]] / h[j]) / 2 / dy[j];
				nmfac[j] = -v1 / 2 / h[j];

				//viscos term
				visnp = dy2h1[j];
				visnc = dy2h2[j];
				visnm = dy2h3[j];

				npfac[j] -= visnp / re;
				ncfac[j] -= visnc / re;
				nmfac[j] -= visnm / re;

				npfac[j] = npfac[j] / 2 * dt;
				ncfac[j] = ncfac[j] / 2 * dt + 1;
				nmfac[j] = nmfac[j] / 2 * dt;
			}
			tdma(nmfac, ncfac, npfac, &dw(i, 1, k), &r3(i, 1, k), n3, n3, n2);
		}

	//solve in x direction
	for (size_t j = 1; j < n2 + 1; j++)
		for (size_t k = 1; k < n3 + 1; k++)
		{
			for (size_t i = 1; i < n1 + 1; i++)
			{
				//d_u^n_w^np_dx
				u1 = (u(i, j, kp[k]) + u(i, j, k)) / 2;
				u2 = (u(ip[i], j, kp[k]) + u(ip[i], j, k)) / 2;
				npfac[i] = u2 / 2 / dx;
				ncfac[i] = (u2 - u1) / 2 / dx;
				nmfac[i] = -u1 / 2 / dx;

				//viscos term
				npfac[i] -= 1 / dx / dx / re;
				ncfac[i] += 2 / dx / dx / re;
				nmfac[i] -= 1 / dx / dx / re;

				npfac[i] = npfac[i] / 2 * dt;
				ncfac[i] = ncfac[i] / 2 * dt + 1;
				nmfac[i] = nmfac[i] / 2 * dt;
			}
			ctdma(nmfac, ncfac, npfac, &dw(1, j, k), &dw(1, j, k), (n2 + 2) * n3, (n2 + 2) * n3, n1);
		}

	//solve in z direction
	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t j = 1; j < n2 + 1; j++)
		{
			for (size_t k = 1; k < n3 + 1; k++)
			{
				//d_w^n_w^np_dz
				w1 = (w(i, j, km[k]) + w(i, j, k)) / 2;
				w2 = (w(i, j, kp[k]) + w(i, j, k)) / 2;
				npfac[k] = w2 / 2 / dz * dt;
				ncfac[k] = (w2 - w1) / 2 / dz * dt;
				nmfac[k] = -w1 / 2 / dz * dt;

				//viscos term
				npfac[k] -= 1 / dz / dz / 2 / re * dt;
				ncfac[k] += 1 / dz / dz / re * dt;
				nmfac[k] -= 1 / dz / dz / 2 / re * dt;

				ncfac[k] += 1;
			}
			ctdma(nmfac, ncfac, npfac, &dw(i, j, 1), &dw(i, j, 1), 1, 1, n3);
		}
}

void form_rp()
{
	getdiv(du, dv, dw, rp, t + dt);
	rp /= dt;
}

void solve_dp()
{
	double alpha, beta;
	double fac;
	int n, m;

	dp.setZero();
	dps_im.setZero();
	dps_re.setZero();
	fft(rp, divs_re, divs_im);
	alpha = 2 * PI / lx;
	beta = 2 * PI / lz;
	for (size_t i = 1; i < n1 + 1; i++)
	{
		m = i - 1 - n1 / 2;
		for (size_t k = 1; k < n3 + 1; k++)
		{
			n = k - 1 - n3 / 2;
			for (size_t j = 1; j < n2 + 1; j++)
			{
				fac = 2 * (cos(alpha * m * dx) - 1) / dx / dx + 2 * (cos(beta * n * dz) - 1) / dz / dz;

				ppfac_re[j] = dy2h1[j];
				pcfac_re[j] = dy2h2[j] + fac;
				pmfac_re[j] = dy2h3[j];
			}
			pcfac_re[n2] = -1 / h[n2] / dy[n2] + fac;
			pmfac_re[n2] = 1 / h[n2] / dy[n2];
			pcfac_re[1] = -1 / h[2] / dy[1] + fac;
			ppfac_re[1] = 1 / h[2] / dy[1];

			pmfac_re[1] = 0;
			ppfac_re[n2] = 0;

			memset(ppfac_im, 0, sizeof(double) * (n2 + 1));
			memset(pcfac_im, 0, sizeof(double) * (n2 + 1));
			memset(pmfac_im, 0, sizeof(double) * (n2 + 1));

			if (m == 0 && n == 0)
			{
				pcfac_re[1] = 1;
				ppfac_re[1] = 0;
				divs_re(i, 1, k) = 0;
				divs_im(i, 1, k) = 0;
			}
			tdma(pmfac_re, pcfac_re, ppfac_re, &dps_re(i, 1, k), &dps_im(i, 1, k), &divs_re(i, 1, k), &divs_im(i, 1, k), n3, n3, n2);
		}
	}

	ifft(dps_re, dps_im, dp);
}

void update_up()
{
	double dpgx, dpgz;
	double dpdy;

	check_flow_rate();

	//update u velocity
	//keep constant mass flow rate
	dpgx = 0;
	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t j = 1; j < n2 + 1; j++)
			for (size_t k = 1; k < n3 + 1; k++)
				dpgx += (du(i, j, k) * dy[j] * dx * dz - dt * (dp(i, j, k) - dp(im[i], j, k)) * dy[j] * dz);

	dpgx = (dpgx - xflow * lx) / lx / ly / lz * dt;

	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t k = 1; k < n3 + 1; k++)
		{
			for (size_t j = 1; j < n2 + 1; j++)
				u(i, j, k) = du(i, j, k) - dt * dpgx - dt * (dp(i, j, k) - dp(im[i], j, k)) / dx;
			u(i, 0, k) = du(i, 0, k);
			u(i, n2 + 1, k) = du(i, n2 + 1, k);
		}

	//update v velocity
	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t j = 2; j < n2 + 1; j++)
			for (size_t k = 1; k < n3 + 1; k++)
			{
				dpdy = (dp(i, j, k) - dp(i, jm[j], k)) / h[j];
				v(i, j, k) = dv(i, j, k) - dt * dpdy;
			}

	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t k = 1; k < n3 + 1; k++)
		{
			v(i, 1, k) = dv(i, 1, k);
			v(i, n2 + 1, k) = dv(i, n2 + 1, k);
		}

	//update w velocity
	//keep 0 flow rate
	dpgz = 0;
	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t j = 1; j < n2 + 1; j++)
			for (size_t k = 1; k < n3 + 1; k++)
			{
				dpgz += (dw(i, j, k) * dy[j] * dx * dz - dt * (dp(i, j, k) - dp(i, j, km[k])) * dx * dy[j]);
			}

	dpgz = (dpgz - zflow * lz) / lx / ly / lz / dt;

	for (size_t i = 1; i < n1 + 1; i++)
		for (size_t k = 1; k < n3 + 1; k++)
		{
			for (size_t j = 1; j < n2 + 1; j++)
			{
				w(i, j, k) = dw(i, j, k) - dt*dpgz - dt * (dp(i, j, k) - dp(i, j, km[k])) / dz;
			}
			w(i, 0, k) = dw(i, 0, k);
			w(i, n2 + 1, k) = dw(i, n2 + 1, k);
		}

	//update pressure
	p += dp;
	pgx += dpgx;
	pgz += dpgz;
}

void output()
{}