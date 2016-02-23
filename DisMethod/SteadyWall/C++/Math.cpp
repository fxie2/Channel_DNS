#include "Math.h"
#include "Datatype.h"
#include "mkl_dfti.h"
#include <iostream>
#include <memory>

//supplimental vector
double *cc, *bb, *dd;
double *cc_im, *bb_im, *dd_im;
double *uu, *vv;
int length = -1;
bool all_allocated = false;

//fft/ifft supplimental array
double *pf_1d_re, *pf_1d_im;
double *sf_1d_re, *sf_1d_im;
Field3d<double> sf2_re, sf2_im;
Field3d<double> shift_re, shift_im;
int fft_dim[4], fft_len = -1;
DFTI_DESCRIPTOR_HANDLE handler;
bool fft_allocated = false;

void initvec(int n)
{
	if (all_allocated)
	{
		if (n <= length) goto init;
		delete[] bb, cc, dd, uu, vv;
		delete[] bb_im, cc_im, dd_im;
		bb = new double[n + 1];
		cc = new double[n + 1];
		dd = new double[n + 1];
		bb_im = new double[n + 1];
		cc_im = new double[n + 1];
		dd_im = new double[n + 1];
		uu = new double[n + 1];
		vv = new double[n + 1];
		length = n;
	}
	else
	{
		bb = new double[n + 1];
		cc = new double[n + 1];
		dd = new double[n + 1];
		bb_im = new double[n + 1];
		cc_im = new double[n + 1];
		dd_im = new double[n + 1];
		uu = new double[n + 1];
		vv = new double[n + 1];
		length = n;
		all_allocated = true;
	}
init:
	memset(bb, 0, sizeof(double) * (n + 1));
	memset(cc, 0, sizeof(double) * (n + 1));
	memset(dd, 0, sizeof(double) * (n + 1));
	memset(bb_im, 0, sizeof(double) * (n + 1));
	memset(cc_im, 0, sizeof(double) * (n + 1));
	memset(dd_im, 0, sizeof(double) * (n + 1));
	memset(uu, 0, sizeof(double) * (n + 1));
	memset(vv, 0, sizeof(double) * (n + 1));
}

void deletevec()
{
	if (all_allocated)
	{
		delete[] bb, cc, dd, uu, vv;
		delete[] bb_im, cc_im, dd_im;
		length = 0;
	}
}

void tdma(double *a, double * b, double * c, double *x, double *r, int gap1, int gap2, int n)
{
	//a, b, c start from index-1
	//x, r start from index-0 with gap1 & gap2
	initvec(n);

	cc[1] = c[1] / b[1];
	for (size_t i = 2; i < n; i++)
	{
		cc[i] = c[i] / (b[i] - a[i] * cc[i - 1]);
	}

	dd[1] = r[0] / b[1];
	for (size_t i = 2; i < n+1; i++)
	{
		dd[i] = (r[(i - 1) * gap2] - a[i] * dd[i - 1]) / (b[i] - a[i] * cc[i - 1]);
	}

	x[(n - 1) * gap1] = dd[n];
	for (size_t i = n - 1; i > 0; i--)
	{
		x[(i - 1) * gap1] = dd[i] - cc[i] * x[i * gap1];
	}
}

void tdma(double *a, double* b, double *c, double *x_re, double *x_im, double *r_re, double *r_im, int gap1, int gap2, int n)
{
	//initvec(n);

	////suppose a, b, c are real
	//cc[1] = c[1] / b[1];
	//for (size_t i = 2; i < n; i++)
	//{
	//	cc[i] = c[i] / (b[i] - a[i] * cc[i - 1]);
	//}

	//dd[1] = r_re[0] / b[1];
	//dd_im[1] = r_im[0] / b[1];
	//for (size_t i = 2; i < n+1; i++)
	//{
	//	dd[i] = (r_re[(i - 1) * gap2] - a[i] * cc[i - 1]) / (b[i] - a[i] * cc[i - 1]);
	//	dd_im[i] = (r_im[(i - 1) * gap2] - a[i] * cc[i - 1]) / (b[i] - a[i] * cc[i - 1]);
	//}

	//x_re[(n - 1) * gap1] = dd[n];
	//x_im[(n - 1) * gap1] = dd_im[n];
	//for (size_t i = n - 1; i > 0; i--)
	//{
	//	x_re[(i - 1) * gap1] = dd[i] - cc[i] * x_re[i * gap1];
	//	x_im[(i - 1) * gap1] = dd_im[i] - cc[i] * x_im[i * gap1];
	//}
	tdma(a, b, c, x_re, r_re, gap1, gap2, n);
	tdma(a, b, c, x_im, r_im, gap1, gap2, n);
}

void ctdma(double *a, double *b, double *c, double *x, double *r, int gap1, int gap2, int n)
{
	initvec(n);

	cc[1] = c[1] / b[1];
	for (size_t i = 2; i < n - 1; i++)
		cc[i] = c[i] / (b[i] - a[i] * cc[i - 1]);

	dd[1] = -a[1] / b[1];
	for (size_t i = 2; i < n - 1; i++)
		dd[i] = -a[i] * dd[i - 1] / (b[i] - a[i] * cc[i - 1]);
	dd[n - 1] = (-c[n - 1] - a[n - 1] * dd[n - 2]) / (b[n - 1] - a[n - 1] * cc[n - 2]);

	uu[n - 1] = dd[n - 1];
	for (size_t i = n - 2; i > 0; i--)
		uu[i] = dd[i] - cc[i] * uu[i + 1];

	dd[1] = r[0] / b[1];
	for (size_t i = 2; i < n; i++)
		dd[i] = (r[(i - 1) * gap2] - a[i] * dd[i - 1]) / (b[i] - a[i] * cc[i - 1]);

	vv[n - 1] = dd[n - 1];
	for (size_t i = n - 2; i > 0; i--)
		vv[i] = dd[i] - cc[i] * vv[i + 1];

	double xn = (r[(n - 1) * gap2] - c[n] * vv[1] - a[n] * vv[n - 1]) /
		(c[n] * uu[1] + a[n] * uu[n - 1] + b[n]);

	for (size_t i = 0; i < n - 1; i++)
		x[i * gap1] = uu[i + 1] * xn + vv[i + 1];
	x[(n - 1) * gap1] = xn;
}

void fft_init(int l1, int l2, int l3)
{
	if (!fft_allocated)
	{
		pf_1d_re = new double[l1 * l3];
		pf_1d_im = new double[l1 * l3];
		sf_1d_re = new double[l1 * l3];
		sf_1d_im = new double[l1 * l3];
		shift_im.reAlloc(l1, l2, l3);
		shift_re.reAlloc(l1, l2, l3);
		sf2_im.reAlloc(l1, l2, l3);
		sf2_re.reAlloc(l1, l2, l3);
		MKL_LONG n[2] = { l1, l3 };
		DftiCreateDescriptor(&handler, DFTI_DOUBLE, DFTI_COMPLEX, 2, n);
		DftiSetValue(handler, DFTI_COMPLEX_STORAGE, DFTI_REAL_REAL);
		DftiSetValue(handler, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
		DftiCommitDescriptor(handler);
		fft_allocated = true;

		fft_dim[1] = l1;
		fft_dim[2] = l2;
		fft_dim[3] = l3;
		fft_len = l1 * l3;
		memset(pf_1d_im, 0, sizeof(double) * l1 * l3);
		memset(pf_1d_re, 0, sizeof(double) * l1 * l3);
		memset(sf_1d_im, 0, sizeof(double) * l1 * l3);
		memset(sf_1d_re, 0, sizeof(double) * l1 * l3);
	}
}

void fft_del()
{
	if (fft_allocated)
	{
		delete[] pf_1d_im, pf_1d_re, sf_1d_im, sf_1d_re;
		DftiFreeDescriptor(&handler);
		fft_allocated = false;
		fft_dim[1] = -1;
		fft_dim[2] = -1;
		fft_dim[3] = -1;
	}
}

void fft(Field3d<double>& pf, Field3d<double>& sf_re, Field3d<double>& sf_im)
{
	fft_init(pf.dim1(), pf.dim2(), pf.dim3());
	for (size_t j = 1; j < pf.dim2() + 1; j++)
	{
		int iter = 0;
		for (size_t i = 1; i <= pf.dim1(); i++)
			for (size_t k = 1; k <= pf.dim3(); k++)
				pf_1d_re[iter++] = pf(i, j, k);
		DftiComputeForward(handler, pf_1d_re, pf_1d_im, sf_1d_re, sf_1d_im);
		int iter2 = 0;
		iter = 0;
		for (size_t i = 1; i <= sf_re.dim1(); i++)
			for (size_t k = 1; k <= sf_re.dim3(); k++)
			{
				sf_re(i, j, k) = sf_1d_re[iter++];
				sf_im(i, j, k) = sf_1d_im[iter2++];
			}
	}
	fftshift(sf_re, sf_im);
}

void ifft(Field3d<double>& sf_re, Field3d<double>& sf_im, Field3d<double>& pf)
{
	fft_init(pf.dim1(), pf.dim2(), pf.dim3());
	sf2_im = sf_im;
	sf2_re = sf_re;
	ifftshift(sf2_re, sf2_im);
	int iter = 0;
	int iter2 = 0;
	for (size_t j = 1; j <= pf.dim2(); j++)
	{
		iter = iter2 = 0;
		for (size_t i = 1; i <= pf.dim1(); i++)
			for (size_t k = 1; k <= pf.dim3(); k++)
			{
				sf_1d_re[iter++] = sf2_re(i, j, k);
				sf_1d_im[iter2++] = sf2_im(i, j, k);
			}
		DftiComputeBackward(handler, sf_1d_re, sf_1d_im, pf_1d_re, pf_1d_im);
		iter = iter2 = 0;
		for (size_t i = 1; i <= pf.dim1(); i++)
			for (size_t k = 1; k <= pf.dim3(); k++)
			{
				pf(i, j, k) = pf_1d_re[iter++];
			}
	}
	pf /= (pf.dim1() * pf.dim3());
}

void fftshift(Field3d<double>& sf_re, Field3d<double>& sf_im)
{
	int n1 = sf_re.dim1();
	int n2 = sf_re.dim2();
	int n3 = sf_re.dim3();
	int n1h = n1 / 2;
	int n3h = n3 / 3;

	for (size_t i = 1; i <= n1h; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3; k++)
				shift_re(i, j, k) = sf_re(n1 - n1h + i, j, k);

	for (size_t i = 1; i <= n1 - n1h; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3; k++)
				shift_re(n1h + i, j, k) = sf_re(i, j, k);

	for (size_t i = 1; i <= n1; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3h; k++)
				sf_re(i, j, k) = shift_re(i, j, n3 - n3h + k);

	for (size_t i = 1; i <= n1; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3 - n3h; k++)
				sf_re(i, j, k + n3h) = shift_re(i, j, k);

	for (size_t i = 1; i <= n1h; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3; k++)
				shift_im(i, j, k) = sf_im(n1 - n1h + i, j, k);

	for (size_t i = 1; i <= n1 - n1h; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3; k++)
				shift_im(n1h + i, j, k) = sf_im(i, j, k);

	for (size_t i = 1; i <= n1; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3h; k++)
				sf_im(i, j, k) = shift_im(i, j, n3 - n3h + k);

	for (size_t i = 1; i <= n1; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3 - n3h; k++)
				sf_im(i, j, k + n3h) = shift_im(i, j, k);

}

void ifftshift(Field3d<double>& sf_re, Field3d<double>& sf_im)
{
	int n1 = sf_re.dim1();
	int n2 = sf_re.dim2();
	int n3 = sf_re.dim3();
	int n1h = n1 / 2;
	int n3h = n3 / 3;

	for (size_t i = 1; i <= n1; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3 - n3h; k++)
				shift_im(i, j, k) = sf_im(i, j, k + n3h);

	for (size_t i = 1; i <= n1; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3h; k++)
				shift_im(i, j, k + n3 - n3h) = sf_im(i, j, k);

	for (size_t i = 1; i <= n1 - n1h; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3; k++)
				sf_im(i, j, k) = shift_im(i + n1h, j, k);

	for (size_t i = 1; i <= n1h; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3; k++)
				sf_im(i + n1 - n1h, j, k) = shift_im(i, j, k);

	for (size_t i = 1; i <= n1; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3 - n3h; k++)
				shift_re(i, j, k) = sf_re(i, j, k + n3h);

	for (size_t i = 1; i <= n1; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3h; k++)
				shift_re(i, j, k + n3 - n3h) = sf_re(i, j, k);

	for (size_t i = 1; i <= n1 - n1h; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3; k++)
				sf_re(i, j, k) = shift_re(i + n1h, j, k);

	for (size_t i = 1; i <= n1h; i++)
		for (size_t j = 1; j <= n2; j++)
			for (size_t k = 1; k <= n3; k++)
				sf_re(i + n1 - n1h, j, k) = shift_re(i, j, k);

}