#include <iostream>
#include <fstream>
#include "Parameter.h"
#include "Field.h"
#include "Mesh.h"
#include "Math.h"

int main()
{
	using namespace std;
	init_parameters();

	//init mesh
	init_mesh();
	//init field
	alloc_field();
	iniup();
	check_div();
	cout << "divmax : " << divmax << '\n';
	//ofstream log;
	//log.open("log.txt", ios::out | ios::trunc);

	while (curnt_step_num < total_step_num)
	{
		solveup();

		if (curnt_step_num % prnt_period == 0)
		{
			check_div();
			check_flow_rate();

			cout << "time :" << t << '\n';
			cout << "divmax : " << divmax << '\n';
			cout << "pgx : " << pgx << '\n';
			cout << "pgz : " << pgz << '\n';
			cout << "xflow : " << xflow << '\n';
			cout << "zflow : " << zflow << '\n';

			cin.get();
			//log << "time :" << t << '\n';
			//log << "divmax : " << divmax << '\n';
			//log << "pgx : " << pgx << '\n';
			//log << "pgz : " << pgz << '\n';
			//log << "xflow : " << xflow << '\n';
			//log << "zflow : " << zflow << '\n';
		}

		curnt_step_num += 1;
		t += dt;
	}

	//log.close();

	del_mesh();
	//dealloc_field();
	deletevec();
	fft_del();

	return 0;
}