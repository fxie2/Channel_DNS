#pragma once
#include "Datatype.h"

typedef Field3d<double> Field;

extern Field u, v, w, p;
extern Field dive;

extern double pgx, pgz;
extern double xflow, zflow;
extern double divmax;
extern double cflmax;

// Allocate field 
void alloc_field();

void dealloc_field();

void iniup();

void check_flow_rate();

void check_div();

void getdiv(Field& uu, Field& vv, Field& ww, Field& ddiv, double tt);

//void check_cfl();

void solveup();

void getvel();

void update_bc();

void getpre();

void getu();

void getv();

void getw();

void finish_vel();

void form_r1();

void form_r2();

void form_r3();

void solve_du();

void solve_dv();

void solve_dw();

void form_rp();

void solve_dp();

void update_up();

void output();