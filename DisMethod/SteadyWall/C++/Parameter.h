#pragma once
#include <string>

using namespace std;

const double PI = 3.141592653589793;

// Geometry parameter
extern size_t n1, n2, n3;
extern double lx, ly, lz;
extern bool use_default_length;
const double DEFAULT_LX = PI;
const double DEFAULT_LY = 2.0;
const double DEFAULT_LZ = PI / 2.0;

// Time parameter
extern double dt;
extern double start_time, end_time;
extern double t;
extern int total_step_num;
extern int curnt_step_num;
extern bool adapted_time_step;

// IO parameter
extern bool init_from_file;
extern bool save_step_data;
extern bool prnt_step_info;
extern string init_file_path;
extern string grid_file_path;
extern string save_file_path;
extern string log_file_path;
extern int prnt_period;
extern int save_period;
extern int insf_save_period;

// Field parameter
extern double re;
extern double init_turb_intensity;
const double DEFAULT_TURB_INTENSITY = 0.01;

// Solve parameter

// Boundary parameter

void init_parameters();

void use_default_parameter();