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
#include "common.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>

static int n_remove;
static gsl_matrix *xmat;

static void compute_xmat(void)
{
  int ii;
  xmat=gsl_matrix_alloc(glob_n_nu,n_remove);
  
  for(ii=0;ii<glob_n_nu;ii++) {
    int jj;
    double nu=glob_nu_arr[ii];
    double lnu=log10(nu);
    double element=1.0;
    for(jj=0;jj<n_remove;jj++) {
      gsl_matrix_set(xmat,ii,jj,element);
      element*=lnu;
    }
  }
}

static void end_polog(void) {
  gsl_matrix_free(xmat);
}

typedef struct {
  double *los;
  double *inv_sigma2;
  double *lnu;
} par_minim;

double model_fit(const gsl_vector *v_c,double lnu)
{
  int jj;
  double polyn=1;
  double exponent=0;

  for(jj=0;jj<n_remove;jj++) {
    double c_j=gsl_vector_get(v_c,jj);
    exponent+=c_j*polyn;
    polyn*=lnu;
  }

  return pow(10.0,exponent);
}

double chi2_func(const gsl_vector *v_c,void *dummy)
{
  int ii;
  par_minim *par=(par_minim *)dummy;
  double chi2=0;

  for(ii=0;ii<glob_n_nu;ii++) {
    double dif;
    double lnu=par->lnu[ii];
    double y=par->los[ii];
    double isig2=par->inv_sigma2[ii];
    double f=model_fit(v_c,lnu);
    dif=f-y;
    chi2+=dif*dif*isig2;
  }
  
  return chi2;
}

static void subtract_foregrounds_nlin(void)
{
#pragma omp parallel default(none)			\
  shared(glob_data_maps,glob_clean_maps,xmat,glob_rms_arr)	\
  shared(glob_nu_arr,n_remove,glob_nth,glob_n_nu)		\
  shared(gsl_multimin_fminimizer_nmsimplex2)
  {
    long ii;
    gsl_multimin_function min_chi2;
    const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *minimizer=gsl_multimin_fminimizer_alloc(T,n_remove);
    gsl_vector *v_c=gsl_vector_alloc(n_remove);
    gsl_vector *step=gsl_vector_alloc(n_remove);
    double *lnu_arr=(double *)my_malloc(glob_n_nu*sizeof(double));
    double *los_minim=(double *)my_malloc(glob_n_nu*sizeof(double));
    double *inv_sigma2_arr=(double *)my_malloc(glob_n_nu*sizeof(double));
    par_minim par;
    
    for(ii=0;ii<glob_n_nu;ii++) lnu_arr[ii]=log10(glob_nu_arr[ii]);

    par.los=los_minim;
    par.lnu=lnu_arr;
    par.inv_sigma2=inv_sigma2_arr;

    min_chi2.n=n_remove;
    min_chi2.f=&chi2_func;
    min_chi2.params=&par;

    for(ii=0;ii<n_remove;ii++) gsl_vector_set(step,ii,0.01);

#pragma omp for
    for(ii=0;ii<glob_nth;ii++) {
      int jj;
      double *los_dirty=&(glob_data_maps[ii*glob_n_nu]);
      double *los_clean=&(glob_clean_maps[ii*glob_n_nu]);

      //Fix the <0 problem
      double tmin=los_dirty[0];
      for(jj=0;jj<glob_n_nu;jj++) {
	double t=los_dirty[jj];
	if(t<tmin) tmin=t;
      }

      double mean_logt=0;
      for(jj=0;jj<glob_n_nu;jj++) {
	double t=los_dirty[jj];
	if(tmin<=0) t+=-2*tmin;
	
	mean_logt+=log10(t);

	inv_sigma2_arr[jj]=1.0/(glob_rms_arr[jj]*glob_rms_arr[jj]);
	los_minim[jj]=t;
      }
      mean_logt/=glob_n_nu;

      gsl_vector_set(v_c,0,mean_logt);
      gsl_vector_set(v_c,1,-2);
      for(jj=2;jj<n_remove;jj++) gsl_vector_set(v_c,jj,0);

      gsl_multimin_fminimizer_set(minimizer,&min_chi2,v_c,step);

      int iter=0;
      int status;
      double size;
      do {
	iter++;
	status=gsl_multimin_fminimizer_iterate(minimizer);
	if(status) break;
	size=gsl_multimin_fminimizer_size(minimizer);
	status=gsl_multimin_test_size(size,1E-6);
      } while((status==GSL_CONTINUE)||(iter<1000));
      
      gsl_vector_memcpy(v_c,minimizer->x);

      for(jj=0;jj<glob_n_nu;jj++) {
	double tfg=model_fit(v_c,lnu_arr[jj]);
	if(tmin<=0) tfg+=2*tmin;
	los_clean[jj]=los_dirty[jj]-tfg;
      }
    }

    gsl_multimin_fminimizer_free(minimizer);
    gsl_vector_free(v_c);
    gsl_vector_free(step);
    free(lnu_arr);
    free(los_minim);
    free(inv_sigma2_arr);
  }
}

static void subtract_foregrounds_lin(void)
{
#pragma omp parallel default(none)			\
  shared(glob_data_maps,glob_clean_maps,xmat,glob_rms_arr)	\
  shared(n_remove,glob_nth,glob_n_nu)
  {
    long ii;
    gsl_multifit_linear_workspace *wsp=
      gsl_multifit_linear_alloc(glob_n_nu,n_remove);
    gsl_vector *v_y=gsl_vector_alloc(glob_n_nu);
    gsl_vector *v_w=gsl_vector_alloc(glob_n_nu);
    gsl_vector *v_c=gsl_vector_alloc(n_remove);
    gsl_matrix *cov_c=gsl_matrix_alloc(n_remove,n_remove);
    gsl_matrix *x_shr=gsl_matrix_alloc(glob_n_nu,n_remove);
    gsl_matrix_memcpy(x_shr,xmat);

#pragma omp for
    for(ii=0;ii<glob_nth;ii++) {
      int jj;
      double chi2;
      double *los_dirty=&(glob_data_maps[ii*glob_n_nu]);
      double *los_clean=&(glob_clean_maps[ii*glob_n_nu]);

      //Fix the <0 problem
      double tmin=los_dirty[0];
      for(jj=0;jj<glob_n_nu;jj++) {
	double t=los_dirty[jj];
	if(t<tmin) tmin=t;
      }

      for(jj=0;jj<glob_n_nu;jj++) {
	double t=los_dirty[jj];

	if(tmin<=0) t+=-2*tmin;

	double inv_sigma=M_LN10*t/glob_rms_arr[jj];
	gsl_vector_set(v_y,jj,log10(t));
	gsl_vector_set(v_w,jj,inv_sigma*inv_sigma);
      }
      gsl_multifit_wlinear(x_shr,v_w,v_y,v_c,cov_c,&chi2,wsp);
      gsl_blas_dgemv(CblasNoTrans,1,x_shr,v_c,0,v_y); //v_y=X*c
      
      for(jj=0;jj<glob_n_nu;jj++) {
	double tfg=pow(10.,gsl_vector_get(v_y,jj));
	if(tmin<=0) tfg+=2*tmin;
	los_clean[jj]=los_dirty[jj]-tfg;
      }
    }

    gsl_multifit_linear_free(wsp);
    gsl_vector_free(v_y);
    gsl_vector_free(v_w);
    gsl_vector_free(v_c);
    gsl_matrix_free(cov_c);
    gsl_matrix_free(x_shr);
  }
}

void do_polog(int nrm)
{
  int use_nlin=0;

  printf("Doing LOS fitting\n");
  n_remove=nrm;
  printf(" - Computing function matrix\n");
  compute_xmat();
  printf(" - Solving least squares and subtracting\n");
  if(use_nlin)
    subtract_foregrounds_nlin();
  else
    subtract_foregrounds_lin();
  end_polog();
  printf("\n");
  postprocess_clean_maps();
}
