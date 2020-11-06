#ifndef _DLL_H_
#define _DLL_H_

#if BUILDING_DLL
#define DLLIMPORT __declspec(dllexport)
#else
#define DLLIMPORT __declspec(dllimport)
#endif

void DLL_EXPORT CTN(short *is_stop, unsigned int* num_monitoring, double* S_data, double* J_data, double* aggregate_number, double* mean_size, double* std_size, double* n_star_data, double* concentration_monomer); //prints out a "Hello World"


#endif
