#ifndef __DLL_H__
#define __DLL_H__

#define DLL_EXPORT extern "C" __declspec(dllexport)
//#define DLL_EXPORT __declspec(dllexport)

DLL_EXPORT void CTN(short *is_stop, unsigned int* num_monitoring, double* S_data, double* J_data, double* aggregate_number, double* mean_size, double* std_size, double* n_star_data, double* monomer_attachment_ns); //prints out a "Hello World"

#endif // __DLL_H__
