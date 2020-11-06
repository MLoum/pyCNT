// dllmain.cpp : Defines the entry point for the DLL application.

//#define _CRT_SECURE_NO_DEPRECATE
#include "stdafx.h"

#include "CTN.h"
#include <windows.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define M_PI 3.14159265358979323846

double get_data_from_ini_file(char* file_buffer, const char *keyword)
{

	char *pt_substring, *pt_strcpy;
	double value;
	short i, n;
	pt_substring = strstr(file_buffer, keyword);

	//offset the pointer with the length of the keyword + 1 for the '=' character
	pt_substring += strlen(keyword) + 1;
	// search for end of line ('\n') character
	n = 0;
	while (pt_substring[n] != '\n')
		n++;
	//printf("n : %d", n);
	pt_strcpy = (char*)malloc(sizeof(char) * (n + 1));

	for (i = 0; i < n; i++)
	{
		pt_strcpy[i] = pt_substring[i];
	}

	pt_strcpy[n] = '\0';

	value = atof(pt_strcpy);

	free(pt_strcpy);
	return value;
}

typedef struct
{
	double number;
	double size;
}aggregates_type;


void DLL_EXPORT CTN(short *is_stop, unsigned int* num_monitoring, double* S_data, double* J_data, double* aggregate_number, double* mean_size, double* std_size, double* n_star_data, double* concentration_monomer)
{
	/*
	is_stop :
	num_monitoring :
	S_data :
	J_data :
	aggregate_number :
	mean_size :
	std_size :
	n_star_data :
	concentration_monomer :
	*/

	// Read input data from ini file
	long lSize;
	char *file_buffer;
	FILE* fp;

	// constant
	double Na = 6.022E23;
	double kB = 1.38E-23;
	
	// Read data from ini file
	// Copy file content in memory
	fp = fopen("CTN.ini", "r");
	//fp = fopen_s("CTN.ini", "r");
	fseek(fp, 0L, SEEK_END);
	lSize = ftell(fp);
	rewind(fp);

	/* allocate memory for entire content */
	file_buffer = (char *) calloc(1, lSize + 1);

	/* copy the file into the buffer */
	fread(file_buffer, lSize, 1, fp);
	fclose(fp);

	// Internal variable i.e. not obtained from ini file
	unsigned int num_iteration, i;
	double sum, mean, variance, std_deviation;	
	// supersaturation
	double S, log_S, S_minus_1;
	// Critical nucleus size (in nb of molecule)
	double n_star;
	// Nucleation rate in per V_solvent_L per time_step_s
	double J;
	// related to growth rate
	double delta_n, constant_n;


	// get data from INI file
	double time_step_ns = get_data_from_ini_file(file_buffer, "time_step_ns");
	double monitoring_time_step_micros = get_data_from_ini_file(file_buffer, "monitoring_time_step_micros");
	int ratio_time_step_monitoring = (int)(monitoring_time_step_micros / time_step_ns * 1E3);
	double end_time_calculation_ms = get_data_from_ini_file(file_buffer, "end_time_calculation_ms");

	unsigned int nb_of_time_step = (unsigned int)(end_time_calculation_ms / time_step_ns * 1E6);
	unsigned int nb_monitoring_point = (unsigned int)(end_time_calculation_ms / monitoring_time_step_micros * 1E3);

	// TODO comment this value ?
	double V_solvent_L = get_data_from_ini_file(file_buffer, "V_solvent_L");
	// monomer initial concentration in mol/L
	double Cini = get_data_from_ini_file(file_buffer, "Cini");

	double from_particle_to_M = 1 / (Na * V_solvent_L);


	double n_monomer_ini = Cini * V_solvent_L * Na;
	double n_monomer = n_monomer_ini;

	// packing density
	double p = get_data_from_ini_file(file_buffer, "p");
	// v : volume of the monomer (not the aggregate)
	double R_nm = get_data_from_ini_file(file_buffer, "R_nm");  // radius in nm (!= R_hydro)
	double v_monomer_nm3 = 4.0 / 3 * M_PI * pow(R_nm, 3);


	// surface tension in N/m or J/m² (monomer/solvant) 19.5 (±0.6), 18.5 (±0.6) and 17.0 (±0.6) cal mol-1  Å-2 in the solvent condition with a DMSO mole fraction of 0.21, 0.26 and 0.32
	double gamma_sl = get_data_from_ini_file(file_buffer, "gamma_sl");

	// Temperature in Kelvin
	double T = get_data_from_ini_file(file_buffer, "T");
	// Viscosity in Pa.s
	double eta = get_data_from_ini_file(file_buffer, "eta");

	// Hydrodynamic radius in nm
	double R_hydro_nm = get_data_from_ini_file(file_buffer, "R_hydro_nm");

	// Diffusion coefficient in m².s-1    
	double D = kB * T / (6 * M_PI*eta*R_hydro_nm*1E-9);
	//printf("D : %.10e\n", D);


	// monomer equilibrium concentration (solubility) in -> m-3 <-
	double Ce_M = get_data_from_ini_file(file_buffer, "Ce_M");
	double Ce_m3 = Ce_M / from_particle_to_M;
	double is_use_Turnbull = (int)(get_data_from_ini_file(file_buffer, "is_use_Turnbull"));
	// beta for spherical aggregate is 0.514
	double beta = 0.514;
	double lambda_heat_over_kT = pow(v_monomer_nm3, 2.0 / 3.0) * gamma_sl / beta / (kB * T);
	if (is_use_Turnbull)
	{
		// in m-3
		Ce_m3 = 1 / v_monomer_nm3 * exp(-lambda_heat_over_kT);
		Ce_M = Ce_m3 * 1000 / Na;
	}
	else
	{
		Ce_m3 = Ce_M * Na / 1E3;
	}


	double cubic_root_36_pi = pow(36 * M_PI, 1.0 / 3.0);
	double pow_1_6_36_pi = pow(36 * M_PI, 1.0 / 6.0);
	double kB_T = kB * T;

	//printf("v_monomer/p : %.10e\n", v_monomer/p);
	double v_monomer_over_p = v_monomer_nm3 / p;
	double v_over_p_two_third = pow(v_monomer_over_p, 2.0 / 3.0)*1E-18;
	//printf("v_over_p_two_third : %.10e\n", v_over_p_two_third);

	// Related to surface tension energy of aggregates - no units
	double theta = gamma_sl * cubic_root_36_pi * v_over_p_two_third / kB_T;
	//printf("theta : %.10e\n", theta);
	double theta_cube = theta * theta * theta;

	// For Nucleation Rate
	double exponential_prefactor_constant_part = (V_solvent_L*1E-3) * time_step_ns * 1E-9* pow_1_6_36_pi * D * (Ce_m3) / (pow(v_monomer_nm3*v_monomer_nm3 / (p*p), 1.0 / 3.0) * 1E-18 * sqrt(theta));
	//printf("exponential_prefactor_constant_part : %0.10e\n", exponential_prefactor_constant_part);

	// For delta n    
	double delta_n_constant_term = 4 * M_PI * D * Ce_m3  * time_step_ns * 1E-9;
	double constant_pow_delta_n = 3 * v_monomer_nm3 / (4 * M_PI*p);

	Ce_M = Ce_m3 * 1E3 / Na;
	// Current concentration of monomer
	double C = n_monomer * from_particle_to_M;

	// Memory Allocation
	aggregates_type *aggregates = (aggregates_type*)malloc(sizeof(aggregates_type)*nb_of_time_step);
	double nb_aggregates = 0;

	// MAIN LOOP
	////////////
	for (num_iteration = 0; num_iteration < nb_of_time_step; num_iteration++)
	{

		if (*is_stop)
			break;

		//printf("num_iteration : %d\n", num_iteration);

		// Recalculate supersaturation
		//printf("C : %f\n", C);
		//printf("Ce_M : %f\n", Ce_M);
		S = C / Ce_M;
		//printf("S : %f\n", S);
		log_S = log(S);
		S_minus_1 = S - 1;
		n_star = pow(2 * theta / (3 * log(S)), 3);
		//printf("n_star : %.10e\n", n_star);
		// nucleation
		//////////////
		// J : Nucleationr rate : new number of super-critic nucleus of size n_star + 1
		// From J in m-3 s-1 (per cubic meter per second) to time_step and V_solvent_L
		J = exponential_prefactor_constant_part * S * log_S * exp(-4 * theta_cube / (27 * log_S*log_S));

		/*
		printf("log_S*log_S : %.10e\n", log_S*log_S);
		printf("theta_cube : %.10e\n", theta_cube);
		printf("exponential_prefactor_constant_part : %.10e\n", exponential_prefactor_constant_part);
		printf("S * log_S : %.10e\n", S * log_S);
		printf("exp(-4*theta_cube/(27*log_S*log_S)) : %.10e\n", exp(-4*theta_cube/(27*log_S*log_S)));
		printf("J : %.10e\n", J);
		*/

		//printf("nb_new_aggregate : %d\n", nb_new_aggregate);

		aggregates[num_iteration].number = J;
		// Important : when a nucleus appears in the solution (it breaks the metasable state), its size is n_star + 1
		aggregates[num_iteration].size = n_star + 1;

		//printf("n_monomer : %.10e\n", n_monomer);
		n_monomer -= aggregates[num_iteration].number * aggregates[num_iteration].size;

		//printf("n_monomer : %.10e\n", n_monomer);

		// growth
		/////////
		constant_n = delta_n_constant_term * S_minus_1*1E-9;
		for (i = 0; i < num_iteration; i++)
		{
			delta_n = constant_n *  pow(aggregates[i].size*constant_pow_delta_n, 1.0 / 3.0);
			//printf("delta_n_constant_term : %.10e\n", delta_n_constant_term);
			//printf("delta_n : %.10e\n", delta_n);
			aggregates[i].size += delta_n;
			n_monomer -= delta_n * aggregates[i].size;
		}


		// Update concentration (in mol/L)
		C = n_monomer * from_particle_to_M;

		if (num_iteration%ratio_time_step_monitoring == 0)
		{
			printf("num_monitoring %d\n : ", *num_monitoring);
			S_data[*num_monitoring] = S;
			J_data[*num_monitoring] = J;
			printf("S : %.10e\n", S);
			printf("J : %.10e\n", J);
			nb_aggregates = 0;
			for (i = 0; i < num_iteration; i++)
			{
				nb_aggregates += aggregates[i].number;
				//printf("aggregates[i].number : %.10e\n", aggregates[i].number);
			}
			aggregate_number[*num_monitoring] = nb_aggregates;
			printf("nb_aggregates : %.10e\n", nb_aggregates);


			mean = 0;
			for (i = 0; i < num_iteration; i++)
			{
				mean += aggregates[i].number * aggregates[i].size;
			}
			mean /= nb_aggregates;
			printf("mean : %.10e\n", mean);
			sum = 0;
			for (i = 0; i < num_iteration; i++)
				sum += pow(aggregates[i].number * aggregates[i].size - mean, 2);

			variance = sum / (float)nb_aggregates;
			std_deviation = sqrt(variance);
			printf("std_deviation : %.10e\n", std_deviation);
			mean_size[*num_monitoring] = mean;
			std_size[*num_monitoring] = std_deviation;
			n_star_data[*num_monitoring] = n_star;
			printf("n_star_data : %.10e\n", n_star);
			concentration_monomer[*num_monitoring] = C;
			printf("concentration_monomer : %.10e\n", C);
			*num_monitoring += 1;
		}
	}// end while main loop

	printf("end calculation ! \n");

	free(file_buffer);
	free(aggregates);
}


/*
BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
    switch (ul_reason_for_call)
    {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
        break;
    }
    return TRUE;
}
*/

BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID lpvReserved)
{
	switch (fdwReason)
	{
	case DLL_PROCESS_ATTACH:
	{
		break;
	}
	case DLL_PROCESS_DETACH:
	{
		break;
	}
	case DLL_THREAD_ATTACH:
	{
		break;
	}
	case DLL_THREAD_DETACH:
	{
		break;
	}
	}

	/* Return TRUE on success, FALSE on failure */
	return TRUE;
}
