#include "Datatype.h"
//#include "Field.h"
#include "Math.h"
#include <iostream>

int main()
{
	double b[6] = { 0,2,2,2,2,2 };
	double a[6] = { -1, 0, -1, -1, -1 , -1};
	double c[6] = { -1, -1, -1, -1, -1 , -1};
	double x[5];
	double r[5] = { 1,2,3,4,5 };
	ctdma(a, b, c, x, r, 1, 1, 5);

	for (size_t i = 0; i < 5;i++)
		std::cout << x[i] << ' ';
	return 0;
}