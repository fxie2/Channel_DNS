#include "Parameter.h"
#include "Mesh.h"
#include <fstream>

using namespace std;
// Vertical position
double *x, *y, *z;

// Cell length
double *h, *dy2h1, *dy2h2, *dy2h3, *dyh1, *dyh2, *dyh3;
double *dy, *dy2dy1, *dy2dy2, *dy2dy3, *dydy1, *dydy2, *dydy3;
double dx, dz;

// Rank
int *im, *ip;
int *jm, *jp;
int *km, *kp;

// Mesh scale in y direction
double scale;

void init_mesh()
{
	x = new double[n1 + 1];
	y = new double[n2 + 1];
	z = new double[n3 + 1];
	h = new double[n2 + 2];
	dyh1 = new double[n2 + 1];
	dyh2 = new double[n2 + 1];
	dyh3 = new double[n2 + 1];
	dy2h1 = new double[n2 + 1];
	dy2h2 = new double[n2 + 1];
	dy2h3 = new double[n2 + 1];
	dy = new double[n2 + 2];
	dy2dy1 = new double[n2 + 1];
	dy2dy2 = new double[n2 + 1];
	dy2dy3 = new double[n2 + 1];
	dydy1 = new double[n2 + 1];
	dydy2 = new double[n2 + 1];
	dydy3 = new double[n2 + 1];
	im = new int[n1 + 1];
	ip = new int[n1 + 1];
	jm = new int[n2 + 1];
	jp = new int[n2 + 1];
	km = new int[n3 + 1];
	kp = new int[n3 + 1];

	dx = lx / n1;
	dz = lz / n3;
	x[0] = 0;
	z[0] = 0;
	for (size_t i = 0; i < n1 + 1; i++)
	{
		x[i] = dx * i;
	}
	for (size_t i = 0; i < n3 + 1; i++)
	{
		z[i] = dz * i;
	}

	//ifstream mesh_input(grid_file_path, ios::_Nocreate, 1);
	//for (size_t i = 0; i < n2 + 1; i++)
	//{
	//	mesh_input >> y[i];
	//}
	//mesh_input.close();
	for (size_t i = 0; i < n2 + 1; i++)
		y[i] = i / (double)n2 * 2.;

	scale = 2.0 / (y[n2] - y[0]);
	for (size_t i = 0; i < n1 + 1; i++) x[i] *= scale;
	for (size_t j = 0; j < n2 + 1; j++) y[j] *= scale;
	for (size_t k = 0; k < n3 + 1; k++) z[k] *= scale;
	ly = 2.0;
	lx = lx * scale;
	lz = lz * scale;
	for (size_t j = 0; j < n2 + 1; j++) y[j] -= 1;

	for (size_t j = 1; j < n2 + 1; j++) dy[j] = y[j] - y[j - 1];
	dy[0] = 0;
	dy[n2 + 1] = 0;
	for (size_t j = 1; j < n2 + 2; j++) h[j] = (dy[j - 1] + dy[j]) / 2;
	for (size_t i = 0; i < n1 + 1; i++)
	{
		im[i] = i - 1;
		ip[i] = i + 1;
	}
	for (size_t j = 0; j < n2 + 1; j++)
	{
		jm[j] = j - 1;
		jp[j] = j + 1;
	}
	for (size_t k = 0; k < n3 + 1; k++)
	{
		km[k] = k - 1;
		kp[k] = k + 1;
	}
	im[1] = n1;
	ip[n1] = 1;
	km[1] = n3;
	kp[n3] = 1;
	for (size_t j = 1; j < n2 + 1; j++)
	{
		dy2dy1[j] = 1 / h[j] / dy[j];
		dy2dy2[j] = -(1 / dy[j] + 1 / dy[jm[j]]) / h[j];
		dy2dy3[j] = 1 / h[j] / dy[jm[j]];
		dydy1[j] = dy[j] * dy[jm[j]] / (dy[j] + dy[jm[j]]) / dy[j] / dy[j];
		dydy2[j] = dy[j] * dy[jm[j]] / (dy[j] + dy[jm[j]]) * (1 / dy[jm[j]] / dy[jm[j]] - 1 / dy[j] / dy[j]);
		dydy3[j] = dy[j] * dy[jm[j]] / (dy[j] + dy[jm[j]]) / dy[jm[j]] / dy[jm[j]] * -1.0;
		dy2h1[j] = 1.0 / h[jp[j]] / dy[j];
		dy2h2[j] = -(h[jp[j]] + h[j]) / h[jp[j]] / h[j] / dy[j];
		dy2h3[j] = 1.0 / h[j] / dy[j];
		dyh1[j] = h[j] * h[jp[j]] / (h[j] + h[jp[j]]) / h[jp[j]] / h[jp[j]];
		dyh2[j] = h[j] * h[jp[j]] / (h[j] + h[jp[j]]) * (1 / h[j] / h[j] - 1 / h[jp[j]] / h[jp[j]]);
		dyh3[j] = h[j] * h[jp[j]] / (h[j] + h[jp[j]]) / h[j] / h[j] * -1.0;
	}

	do{
		int j;
		j = 1;
		dy2h1[j] = 2 / (h[j] + h[jp[j]]) / h[jp[j]];
		dy2h2[j] = 2 / (h[j] + h[jp[j]]) * -(1 / h[jp[j]] + 1 / h[j]);
		dy2h3[j] = 2 / (h[j] + h[jp[j]]) / h[j];
		j = n2;
		dy2h1[j] = 2 / (h[j] + h[jp[j]]) / h[jp[j]];
		dy2h2[j] = 2 / (h[j] + h[jp[j]]) * -(1 / h[jp[j]] + 1 / h[j]);
		dy2h3[j] = 2 / (h[j] + h[jp[j]]) / h[j];
	} while (0);
}

void del_mesh()
{
	delete[] x, y, z;
	delete[] h, dyh1, dyh2, dyh3, dy2h1, dy2h2, dy2h3;
	delete[] dy, dydy1, dydy2, dydy3, dy2dy1, dy2dy2, dy2dy3;
	delete[] im, ip, jm, jp, km, kp;
}
