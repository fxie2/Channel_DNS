#include "Datatype.h"
//#include "Field.h"
#include "Math.h"
#include <iostream>

int main()
{
	using namespace std;
	Field3d<double> a(3, 1, 3), a0(3, 1, 3);
	Field3d<double> bre(3, 1, 3), bim(3, 1, 3);
	for (size_t i = 1; i < 4; i++)
		for (size_t j = 1; j < 4; j++)
			a(i, 1, j) = i * j;
	a0 = a;
	fft(a, bre, bim);
	cout << 'a' << '\n';
	cout << a;
	cout << "bre\n";
	cout << bre;
	cout << "bim\n";
	cout << bim;
	ifft(bre, bim, a);
	cout << 'a' << '\n';
	cout << a;
	cout << "bre\n";
	cout << bre;
	cout << "bim\n";
	cout << bim;
	cout << "a0 - a\n" << a0 - a;
	return 0;
}