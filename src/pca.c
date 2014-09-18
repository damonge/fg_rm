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
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_blas.h>

static int n_remove;
static gsl_matrix *covariance;
static gsl_vector *eigenvals;
static gsl_matrix *eigenvecs;
static gsl_matrix *fg_removal;

static void compute_covariance(void)
{
  double *cov=my_calloc(glob_n_nu*glob_n_nu,sizeof(double));
  
  //TODO: a different ordering of the maps might make this faster...
#pragma omp parallel default(none)		\
  shared(glob_data_maps,cov,glob_n_nu,glob_nth)
  {
    long ii;
    double norm=1.0/glob_nth;
    double *cov_thr=my_calloc(glob_n_nu*glob_n_nu,sizeof(double));
    
#pragma omp for
    for(ii=0;ii<glob_nth;ii++) {
      int j1;
      double *los=&(glob_data_maps[ii*glob_n_nu]);
      for(j1=0;j1<glob_n_nu;j1++) {
	int j2;
	for(j2=0;j2<glob_n_nu;j2++) {
	  int index=j2+j1*glob_n_nu;
	  cov_thr[index]+=los[j1]*los[j2];
	}
      }
    }

#pragma omp critical
    {
      int i;
      for(i=0;i<glob_n_nu;i++) {
	int j;
	for(j=0;j<glob_n_nu;j++) {
	  int index=j+i*glob_n_nu;
	  cov[index]+=cov_thr[index]*norm;
	}
      }
    }
    
    free(cov_thr);
  }

  covariance=gsl_matrix_alloc(glob_n_nu,glob_n_nu);

  int i;
  for(i=0;i<glob_n_nu;i++) {
    int j;
    for(j=0;j<glob_n_nu;j++) {
      int index=j+glob_n_nu*i;
      double inv_var=1.0/(glob_rms_arr[j]*glob_rms_arr[i]);
      gsl_matrix_set(covariance,i,j,cov[index]*inv_var);
    }
  }
  free(cov);
}

static void diagonalize_covariance(void)
{
  gsl_vector *vec_dum=gsl_vector_alloc(glob_n_nu);
  gsl_matrix *evec_dum=gsl_matrix_alloc(glob_n_nu,glob_n_nu);
  gsl_vector *eval_dum=gsl_vector_alloc(glob_n_nu);
  eigenvals=gsl_vector_alloc(glob_n_nu);
  eigenvecs=gsl_matrix_alloc(glob_n_nu,glob_n_nu);

  //Diagonalize
  gsl_eigen_symmv_workspace *w=gsl_eigen_symmv_alloc(glob_n_nu);
  gsl_eigen_symmv(covariance,eval_dum,evec_dum,w);
  gsl_eigen_symmv_free(w);

  //Sort eigenvalues
  gsl_permutation *p=gsl_permutation_alloc(glob_n_nu);
  gsl_sort_vector_index(p,eval_dum);
  
  int ii;
  for(ii=0;ii<glob_n_nu;ii++) {
    int inew=gsl_permutation_get(p,ii);
    gsl_vector_set(eigenvals,ii,gsl_vector_get(eval_dum,inew));
    gsl_matrix_get_col(vec_dum,evec_dum,inew);
    gsl_matrix_set_col(eigenvecs,ii,vec_dum);
  }
  gsl_permutation_free(p);
  gsl_vector_free(vec_dum);
  gsl_vector_free(eval_dum);
  gsl_matrix_free(evec_dum);

  FILE *fo;
  char fname[256];
  sprintf(fname,"%s_pca_eigvals.dat",glob_prefix_out);
  fo=my_fopen(fname,"w");
  for(ii=0;ii<glob_n_nu;ii++) {
    double lambda=gsl_vector_get(eigenvals,ii);
    fprintf(fo,"%d %lE\n",ii,lambda);
  }
  fclose(fo);
}

static void compute_fg_removal(void)
{
  int ii;
  gsl_matrix *id=gsl_matrix_alloc(glob_n_nu,glob_n_nu);
  gsl_matrix *dum=gsl_matrix_alloc(glob_n_nu,glob_n_nu);
  fg_removal=gsl_matrix_alloc(glob_n_nu,glob_n_nu);
  gsl_matrix_set_identity(id);
  gsl_matrix_set_zero(dum);
  gsl_matrix_set_zero(fg_removal);
  
  //Eliminating the first n_remove principal eigenvectors
  for(ii=0;ii<n_remove;ii++)
    gsl_matrix_set(id,glob_n_nu-ii-1,glob_n_nu-ii-1,0);

  //Computing FG=U*S*U^T
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,id,eigenvecs,0.0,dum);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,eigenvecs,dum,0.0,fg_removal);
  
  gsl_matrix_free(id);
  gsl_matrix_free(dum);
}

static void end_pca(void) {
  gsl_matrix_free(covariance);
  gsl_matrix_free(eigenvecs);
  gsl_vector_free(eigenvals);
}

static void subtract_foregrounds(void)
{
#pragma omp parallel default(none)		\
  shared(glob_data_maps,glob_clean_maps,fg_removal)	\
  shared(glob_rms_arr,glob_nth,glob_n_nu)
  {
    long ii;
    gsl_vector *v_in=gsl_vector_alloc(glob_n_nu);
    gsl_vector *v_out=gsl_vector_alloc(glob_n_nu);
    gsl_matrix *fgrm=gsl_matrix_alloc(glob_n_nu,glob_n_nu);

    gsl_matrix_memcpy(fgrm,fg_removal);

#pragma omp for
    for(ii=0;ii<glob_nth;ii++) {
      int jj;
      double *los_dirty=&(glob_data_maps[ii*glob_n_nu]);
      double *los_clean=&(glob_clean_maps[ii*glob_n_nu]);
      for(jj=0;jj<glob_n_nu;jj++)
	gsl_vector_set(v_in,jj,los_dirty[jj]/glob_rms_arr[jj]);
      gsl_blas_dgemv(CblasNoTrans,1.0,fg_removal,v_in,0,v_out);
      for(jj=0;jj<glob_n_nu;jj++)
	los_clean[jj]=glob_rms_arr[jj]*gsl_vector_get(v_out,jj);
    }

    gsl_vector_free(v_in);
    gsl_vector_free(v_out);
    gsl_matrix_free(fgrm);
  }
}

void do_pca(int nrm)
{
  printf("Doing PCA\n");
  n_remove=nrm;
  printf(" - Computing covariance\n");
  compute_covariance();
  printf(" - Diagonalizing covariance\n");
  diagonalize_covariance();
  printf(" - Computing foreground removal matrix\n");
  compute_fg_removal();
  printf(" - Removing foreground\n");
  subtract_foregrounds();
  printf("\n");
  end_pca();
  postprocess_clean_maps();
}
