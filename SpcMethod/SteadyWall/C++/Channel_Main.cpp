/*********************************************************
\file			Channel_Main.cpp
\author			Luohao Wang
\version		1.0
\date			2015-9-9
\brief			The main program
********************************************************/

#include <stdio.h>
#include "Parameters.h"
#include "Additional_functions.h"

int main(int argc[], char* argv[])
{
	Global_Parameters parm;
	if (0 != initialize(parm)) return -1;
	return 0;
}