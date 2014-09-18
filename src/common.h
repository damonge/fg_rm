///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of fg_rm.                                       //
//                                                                   //
// fg_rm is free software: you can redistribute it and/or modify it  //
// under the terms of the GNU General Public License as published by //
// the Free Software Foundation, either version 3 of the License, or //
// (at your option) any later version.                               //
//                                                                   //
// fg_rm is distributed in the hope that it will be useful, but      //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with fg_rm.  If not, see <http://www.gnu.org/licenses/>.    //
//                                                                   //
///////////////////////////////////////////////////////////////////////
#ifndef _COMMON_H_
#define _COMMON_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

#define NU_21CM 1420.4
#define OMEGA_M 0.3

extern char glob_prefix_out[256];
extern char glob_prefix_cosmo_in[256];
extern char glob_prefix_in[256];
extern char glob_fname_nuTable[256];
extern char glob_fname_mask[256];
extern char glob_fname_rms[256];
extern char glob_method[32];

extern int glob_analyze_pk;
extern int glob_output_maps;

extern int glob_calculate_rms;
extern int glob_n_foregrounds;
extern int glob_start_radial;
extern int glob_end_radial;
extern int glob_n_bins_radial;

extern double glob_fsky;

extern double *glob_nu_arr;
extern double *glob_nu0_arr;
extern double *glob_nuf_arr;
extern int *glob_map_suffixes;
extern double *glob_rms_arr;

extern int glob_n_nu;
extern long glob_nside;
extern long glob_nth;

extern long *glob_indices_in_mask;
extern float *glob_mask;

extern double *glob_data_maps;
extern double *glob_clean_maps;

// Defined in common.c
void report_error(int level,char *fmt,...);
void *my_malloc(size_t size);
void *my_calloc(size_t nmemb,size_t size);
FILE *my_fopen(const char *path,const char *mode);
void read_data(char *fname_params);
void end_all(void);
void write_cleaned_maps();
void postprocess_clean_maps(void);
  
// Defined in cls_analysis.c 
void do_analysis(int i_start,int i_end,int n_superbins);
  
// Defined in pca.c
void do_pca(int nrm);
  
// Defined in polog.c
void do_polog(int nrm);
  
// Defined in fastica.cpp
void do_fastica(int nrm);

#ifdef __cplusplus
}
#endif //__cplusplus

#endif //_COMMON_H_
