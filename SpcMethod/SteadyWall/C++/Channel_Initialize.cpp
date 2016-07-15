/*********************************************************
\file			Channel_Initialize.cpp
\author			Luohao Wang
\version		1.0
\date			2015-9-9
\brief			Initialize parameters
********************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "Parameters.h"

#ifdef WIN32
int initialize_from_file(Global_Parameters &para, char *path = ".\\init.dat")
#endif
#ifdef linux
int initialize_from_file(Global_Parameters &para, char *path = "./init.dat")
#endif
{
	FILE *fp;
	fp = fopen(path, "rb");
	if (NULL != fp)
	{
		fread(&para, sizeof(Global_Parameters), 1, fp);
		return 0;
	}
	else
	{
		return -1;
	}
}

int initialize_manually(Global_Parameters &para)
{
	printf("Please input the following parameters : \n");
	printf("\nMesh size in X : ");		scanf("%d", &(para.real_size.x));
	printf("\nMesh size in Y : ");		scanf("%d", &(para.real_size.y));
	printf("\nMesh size in Z : ");		scanf("%d", &(para.real_size.z));
	printf("\nReynolds number : ");		scanf("%f", &(para.re));
	printf("\nWave number in X : ");	scanf("%f", &(para.alpha));
	printf("\nWave number in Z : ");	scanf("%f", &(para.beta));
	printf("\nPressure gradient : ");	scanf("%f", &(para.dpdx));
	printf("\nMass flux : ");			scanf("%f", &(para.mass_flux));
	printf("\nTime step size : ");		scanf("%f", &(para.dt));
	printf("\nStart time : ");			scanf("%f", &(para.start_time));
	printf("\nEnd time : ");			scanf("%f", &(para.end_time));

	//Set half_size and trans_size value
	para.half_size.x = para.real_size.x / 2;
	para.half_size.y = para.real_size.y / 2;
	para.half_size.z = para.real_size.z / 2;

	para.trans_size.x = para.real_size.x / 2 * 3;
	para.trans_size.y = para.real_size.y / 2 * 3;
	para.trans_size.z = para.real_size.z / 2 * 3;

	//Set step_num value
	para.step_num =  floor(para.end_time - para.start_time) / para.dt;
	return 0;
}

//!Set parameters' values
int initialize(Global_Parameters &para)
{
	int not_finished = 1;
	char command;
	char path[256] = {'\0'};
	while (0 != not_finished)
	{
		printf("Initialize from file or manually? [F/M/Q] : ");
		scanf(" %c", &command);
		switch (command)
		{
		case 'F':
		case 'f':
			printf("Please input the path for the file, ENTER for deafult : \n");
			getchar();
			if ('\n' != (command = getchar()))
			{		
				ungetc(command, stdin);
				scanf("%s", path);
				not_finished = initialize_from_file(para, path);
			}
			else
				not_finished = initialize_from_file(para);
			break;
		case 'M':
		case 'm':
			not_finished = initialize_manually(para);
			break;
		case 'Q':
		case 'q':
			not_finished = 0;
			return -1;
			break;
		default:
			break;
		}
	}
	return 0;
}