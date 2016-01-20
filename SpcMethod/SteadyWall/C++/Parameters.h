/*********************************************************
 \file			Parameters.h
 \author		Luohao Wang
 \version		1.0
 \date			2015-9-9
 \brief			Declear the parameters in the computation
 ********************************************************/

#pragma once

///Struct to store the mesh size in 3-dimension
struct Mesh_size
{
	///Mesh size in 3-dimension
	unsigned int x, y, z;
};

///Struct to store the global parameters
struct Global_Parameters
{
public:
	Mesh_size real_size;		///<Mesh size of the real problem
	Mesh_size half_size;		///<Mesh size of the half problem
	Mesh_size trans_size;	///<Mesh size for 3/2 rule
	double	 re;			///<Reynolds number
	double	 alpha;			///<Wave number in X-direction
	double	 beta;			///<Wave number in Z-direction
	double	 dpdx;			///<Pressure gradient in X-direction
	double	 mass_flux;		///<Mass flux
	double	 dt;			///<Time space
	double	 start_time;	///<Start time
	double	 end_time;		///<End time
	int		 step_num;		///<Time step numbers
};	///<Global parameter struct