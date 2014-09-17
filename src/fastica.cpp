#include <itpp/itsignal.h>
#include "common.h"

using namespace itpp;
using namespace std;

void do_fastica(int nrm)
{
  long ii;

  printf("Sorting data\n");
  mat X=zeros(glob_n_nu,glob_nth);

  double *mean_arr=(double *)my_calloc(glob_n_nu,sizeof(double));
  for(ii=0;ii<glob_nth;ii++) {
    int jj;
    for(jj=0;jj<glob_n_nu;jj++)
      mean_arr[jj]+=glob_data_maps[jj+ii*glob_n_nu];
  }
  for(ii=0;ii<glob_n_nu;ii++)
    mean_arr[ii]/=glob_nth;
  
  for(ii=0;ii<glob_nth;ii++) {
    int jj;
    for(jj=0;jj<glob_n_nu;jj++)
      X(jj,ii)=glob_data_maps[jj+ii*glob_n_nu]-mean_arr[jj];
  }
  

  printf("Setting ICA\n");
  Fast_ICA ica(X);
  ica.set_nrof_independent_components(nrm);
  //  ica.set_approach(FICA_APPROACH_DEFL);
  //  ica.set_fine_tune(false);
  ica.set_max_num_iterations(100);
  //  ica.set_non_linearity(FICA_NONLIN_GAUSS);
  ica.set_epsilon(0.01);

  printf("Doing FASTICA\n");
  ica.separate();

  printf("Getting matrices\n");
  mat A=ica.get_mixing_matrix();
  mat S=ica.get_independent_components();
  mat FG=A*S;

  /*
    printf("X %d %d\n",X.rows(),X.cols());
    printf("A %d %d\n",A.rows(),A.cols());
    printf("FG %d %d\n",FG.rows(),FG.cols());
  */

  printf("Subtracting foregrounds\n");
  for(ii=0;ii<glob_nth;ii++) {
    int jj;
    for(jj=0;jj<glob_n_nu;jj++) {
      glob_clean_maps[jj+ii*glob_n_nu]=glob_data_maps[jj+ii*glob_n_nu]
	-mean_arr[jj]-FG(jj,ii);
    }
  }
  
  free(mean_arr);
}
