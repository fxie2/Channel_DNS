#pragma once
#include "Datatype.h"
//#include "mkl_dfti.h"

void initvec(int n);

void deletevec();

//bugs in start loc of a,b,c, fix it
void tdma(double *a, double * b, double * c, double *x, double *r, int gap1, int gap2, int n);

void tdma(double *a, double* b, double *c, double *x_re, double *x_im, double *r_re, double *r_im, int gap1, int gap2, int n);

void ctdma(double *a, double *b, double *c, double *x, double *r, int gap1, int gap2, int n);

void fft_init(int l1, int l2, int l3);

void fft_del();

void fft(Field3d<double>& pf, Field3d<double>& sf_re, Field3d<double>& sf_im);

void ifft(Field3d<double>& sf_re, Field3d<double>& sf_im, Field3d<double>& pf);

void fftshift(Field3d<double>& sf_re, Field3d<double>& sf_im);

void ifftshift(Field3d<double>& sf_re, Field3d<double>& sf_im);