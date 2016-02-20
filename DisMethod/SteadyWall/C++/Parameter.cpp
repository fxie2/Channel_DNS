#include "Parameter.h"

size_t n1, n2, n3;
double lx, ly, lz;
bool use_default_length;

// Time parameter
double dt;
double start_time, end_time;
double t;
int total_step_num;
int curnt_step_num;
bool adapted_time_step;

// IO parameter
bool init_from_file;
bool save_step_data;
bool prnt_step_info;
string init_file_path;
string grid_file_path;
string save_file_path;
string log_file_path;
int prnt_period;
int save_period;
int insf_save_period;

// Field parameter
double re;
double init_turb_intensity;

void init_parameters()
{
	use_default_parameter();
}

void use_default_parameter()
{
	n1 = 64;
	n2 = 64;
	n3 = 64;
	lx = DEFAULT_LX;
	ly = DEFAULT_LY;
	lz = DEFAULT_LZ;

	dt = 0.01;
	start_time = 0;
	total_step_num = 1000;
	end_time = start_time + dt * total_step_num;
	curnt_step_num = 0;
	t = dt * curnt_step_num;
	adapted_time_step = true;

	init_from_file = false;
	save_step_data = true;
	prnt_step_info = true;
	grid_file_path = "channel.grd";
	save_file_path = "f:\\dataxxxx\\";
	log_file_path = save_file_path;
	prnt_period = 1;
	save_period = 100;
	insf_save_period = 10;

	re = 1000.0;
	init_turb_intensity = DEFAULT_TURB_INTENSITY;
}
