#include "common.h"
#ifdef __cplusplus
extern "C" {
#endif //__cplusplus
#include "healpix_extra.h"
#ifdef __cplusplus
}
#endif //__cplusplus
#include "fftw3.h"

static double *cosmo_maps,*res_maps;

static void end_analysis(void)
{
  free(cosmo_maps);
  free(res_maps);
}

static void init_analysis(void)
{
  int ii;

  cosmo_maps=(double *)my_malloc(glob_n_nu*glob_nth*sizeof(double));
  res_maps=(double *)my_malloc(glob_n_nu*glob_nth*sizeof(double));

  //Read signal map and subtract mean
  for(ii=0;ii<glob_n_nu;ii++) {
    long jj,ns;
    char fname[256];
    float *map_signal;

#ifdef _VERBOSE
    printf("%d \n",ii);
#endif //_VERBOSE

    double mean=0;
    sprintf(fname,"%s_%03d.fits",glob_prefix_cosmo_in,glob_map_suffixes[ii]);
    map_signal=he_read_healpix_map(fname,&ns,0);
    for(jj=0;jj<glob_nth;jj++) {
      double x=(double)(map_signal[glob_indices_in_mask[jj]]);

      mean+=x;
    }
    mean/=glob_nth;
    for(jj=0;jj<12*glob_nside*glob_nside;jj++) {
      if(glob_mask[jj]<0.01)
	map_signal[jj]=0;
      else
	map_signal[jj]-=mean;
    }

    for(jj=0;jj<glob_nth;jj++) {
      double cosmo=map_signal[glob_indices_in_mask[jj]];
      double clean=glob_clean_maps[ii+jj*glob_n_nu];
      cosmo_maps[ii+jj*glob_n_nu]=cosmo;
      res_maps[ii+jj*glob_n_nu]=clean-cosmo;
    }
    free(map_signal);
  }
}

static double *compute_pkrad(double *maps,int i_x0,int n_x)
{
  int i;
  double *pk=(double *)my_calloc((n_x/2+1),sizeof(double));
  
  int nthr=omp_get_max_threads();
  double **loss=(double **)my_malloc(nthr*sizeof(double *));
  fftw_complex **dfts=(fftw_complex **)my_malloc(nthr*sizeof(fftw_complex *));
  fftw_plan *plans=(fftw_plan *)my_malloc(nthr*sizeof(fftw_plan));
  for(i=0;i<nthr;i++) {
    loss[i]=(double *)my_malloc(n_x*sizeof(double));
    dfts[i]=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(n_x/2+1));
    plans[i]=fftw_plan_dft_r2c_1d(n_x,loss[i],dfts[i],FFTW_ESTIMATE);
  }

#pragma omp parallel default(none)			\
  shared(n_x,glob_nth,i_x0,glob_n_nu,pk,maps,dfts,loss,plans)
  {
    long ii;
    int pid=omp_get_thread_num();
    double *pk_thr=(double *)my_calloc((n_x/2+1),sizeof(double));
    
#pragma omp for
    for(ii=0;ii<glob_nth;ii++) {
      int jj;
      for(jj=i_x0;jj<i_x0+n_x;jj++)
      	loss[pid][jj-i_x0]=maps[jj+ii*glob_n_nu];
      fftw_execute(plans[pid]);
      for(jj=0;jj<n_x/2+1;jj++)
	pk_thr[jj]+=creal(dfts[pid][jj]*conj(dfts[pid][jj]));
    }// end omp for

#pragma omp critical
    {
      for(ii=0;ii<n_x/2+1;ii++)
	pk[ii]+=pk_thr[ii]/glob_nth;
    }// end omp critical

    free(pk_thr);
  }// end omp parallel

  for(i=0;i<nthr;i++) {
    free(loss[i]);
    fftw_free(dfts[i]);
    fftw_destroy_plan(plans[i]);
  }
  free(loss);
  free(dfts);
  free(plans);

  return pk;
}

static void analyze_pkrad(int i_start,int i_end,int n_superbins)
{
  int ii;

  //Iterate over frequency "superbins"
  int n_x=(i_end-i_start+1)/n_superbins;
  for(ii=0;ii<n_superbins;ii++) {
#ifdef _VERBOSE
    printf("%d \n",ii);
#endif //_VERBOSE

    FILE *fo;
    int jj;
    char fname[256];
    int i_x0=i_start+ii*n_x;

    double *fft_k_cosmo=compute_pkrad(cosmo_maps,i_x0,n_x);
    double *fft_k_clean=compute_pkrad(glob_clean_maps,i_x0,n_x);
    double *fft_k_res=compute_pkrad(res_maps,i_x0,n_x);
    
    sprintf(fname,"%s_pkrad_%03d.dat",glob_prefix_out,ii);
    fo=my_fopen(fname,"w");
    
    double zmean=NU_21CM/(0.5*(glob_nuf_arr[i_x0+n_x-1]+glob_nu0_arr[i_x0]))-1;
    double nu_interval=glob_nuf_arr[i_x0+n_x-1]-glob_nu0_arr[i_x0];
    double delta_k_nu=2*M_PI/nu_interval;
    double q=0.47366*sqrt(OMEGA_M*(1+zmean)*(1+zmean)*(1+zmean)+(1-OMEGA_M))/
      ((1+zmean)*(1+zmean));
    for(jj=0;jj<n_x/2+1;jj++) {
      double k_r=jj*delta_k_nu*q;
      double pk_cosmo=fft_k_cosmo[jj]*nu_interval/(q*n_x*n_x);
      double pk_res=fft_k_res[jj]*nu_interval/(q*n_x*n_x);
      double pk_clean=fft_k_clean[jj]*nu_interval/(q*n_x*n_x);

      fprintf(fo,"%lE %lE %lE %lE\n",k_r,pk_clean,pk_cosmo,pk_res);
    }
    fclose(fo);

    free(fft_k_cosmo);
    free(fft_k_res);
    free(fft_k_clean);
  }
}

static void analyze_cls()
{
  int ii;
  int lmax=3*glob_nside-1;
  double norm_fsky=1./glob_fsky;
  float *map_cosmo=(float *)my_calloc(12*glob_nside*glob_nside,sizeof(float));
  float *map_clean=(float *)my_calloc(12*glob_nside*glob_nside,sizeof(float));
  float *map_res=(float *)my_calloc(12*glob_nside*glob_nside,sizeof(float));
  float *cls_cosmo=(float *)my_malloc((lmax+1)*sizeof(float));
  float *cls_clean=(float *)my_malloc((lmax+1)*sizeof(float));
  float *cls_res=(float *)my_malloc((lmax+1)*sizeof(float));
 
  for(ii=0;ii<glob_n_nu;ii++) {
    long jj;
    char fname[256];

#ifdef _VERBOSE
    printf("%d \n",ii);
#endif //_VERBOSE

    //Form maps
    for(jj=0;jj<glob_nth;jj++) {
      map_cosmo[glob_indices_in_mask[jj]]=(float)(cosmo_maps[ii+jj*glob_n_nu]);
      map_clean[glob_indices_in_mask[jj]]=(float)(glob_clean_maps[ii+jj*glob_n_nu]);
      map_res[glob_indices_in_mask[jj]]=(float)(res_maps[ii+jj*glob_n_nu]);
    }

    //Compute cls
    he_anafast(glob_nside,lmax,map_cosmo,cls_cosmo);
    he_anafast(glob_nside,lmax,map_clean,cls_clean);
    he_anafast(glob_nside,lmax,map_res,cls_res);

    FILE *fo;
    sprintf(fname,"%s_cls_%03d.dat",glob_prefix_out,glob_map_suffixes[ii]);
    fo=my_fopen(fname,"w");
    for(jj=0;jj<=lmax;jj++) {
      double cl_cosmo=cls_cosmo[jj]*norm_fsky;
      double cl_clean=cls_clean[jj]*norm_fsky;
      double cl_res=cls_res[jj]*norm_fsky;
      double cvar=cl_cosmo*sqrt(2*norm_fsky/(2*jj+1.0));
      fprintf(fo,"%ld %lE %lE %lE %lE\n",jj,
	      cl_clean,cl_cosmo,cl_res,cvar);
    }
    fclose(fo);
  }

  free(map_cosmo);
  free(map_clean);
  free(map_res);
  free(cls_cosmo);
  free(cls_clean);
  free(cls_res);
}

void do_analysis(int i_start,int i_end,int n_superbins)
{
  printf("Reading cosmo map and computing residuals\n");
  init_analysis();
  printf("Analyzing cls\n");
  analyze_cls();
  printf("Analyzing radial pk\n");
  analyze_pkrad(i_start,i_end,n_superbins);
  end_analysis();
}
