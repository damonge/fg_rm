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
#ifdef __cplusplus
extern "C" {
#endif //__cplusplus
#include "healpix_extra.h"
#ifdef __cplusplus
}
#endif //__cplusplus

char glob_prefix_out[256]="output";
char glob_prefix_cosmo_in[256]="signal";
char glob_prefix_in[256]="map";
char glob_fname_nuTable[256]="nutable.txt";
char glob_fname_mask[256]="mask.fits";
char glob_fname_rms[256]="rms.txt";
char glob_method[32]="pca";

int glob_analyze_pk=0;
int glob_output_maps=0;

int glob_calculate_rms=0;
int glob_n_foregrounds=5;
int glob_start_radial=0;
int glob_end_radial=100;
int glob_n_bins_radial=2;

double glob_fsky;

double *glob_nu_arr;
double *glob_nu0_arr;
double *glob_nuf_arr;
int *glob_map_suffixes;
double *glob_rms_arr;

int glob_n_nu;
long glob_nside;
long glob_nth;

long *glob_indices_in_mask;
float *glob_mask;

double *glob_data_maps;
double *glob_clean_maps;

void report_error(int level,char *fmt,...)
{
  va_list args;
  char head[256];
  char msg[256];

  sprintf(head,"Error found in file %s, line %d\n",__FILE__,__LINE__);
  va_start(args,fmt);
  vsprintf(msg,fmt,args);
  va_end(args);

  if(level) {
    fprintf(stderr,"   Fatal: %s   %s",head,msg);
    exit(level);
  }
  else
    fprintf(stderr,"   Warning: %s   %s",head,msg);
}

FILE *my_fopen(const char *path,const char *mode)
{
  FILE *fout=fopen(path,mode);
  if(fout==NULL) 
    report_error(1,"Couldn't open file %s\n",path);

  return fout;
}

void *my_malloc(size_t size)
{
  void *outptr=malloc(size);
  if(outptr==NULL) report_error(1,"Out of memory\n");

  return outptr;
}

void *my_calloc(size_t nmemb,size_t size)
{
  void *outptr=calloc(nmemb,size);
  if(outptr==NULL) 
    report_error(1,"Out of memory\n");

  return outptr;
}

static void error_read_line(char *fname,int line)
{
  report_error(1,"Error reading %s, line %d\n",fname,line);
}

static int linecount(FILE *f)
{
  int i0=0;
  char ch[1024];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}

static void read_nu_arr(void)
{
  int ii;
  FILE *fil=my_fopen(glob_fname_nuTable,"r");
  glob_n_nu=linecount(fil);
  rewind(fil);

  glob_nu_arr=(double *)my_malloc(glob_n_nu*sizeof(double));
  glob_nu0_arr=(double *)my_malloc(glob_n_nu*sizeof(double));
  glob_nuf_arr=(double *)my_malloc(glob_n_nu*sizeof(double));
  glob_map_suffixes=(int *)my_malloc(glob_n_nu*sizeof(int));

  for(ii=0;ii<glob_n_nu;ii++) {
    int suf;
    double nu0,nu1,z0,z1;
    int st=fscanf(fil,"%d %lf %lf %lf %lf\n",&suf,
		  &nu0,&nu1,&z0,&z1);
    if(st!=5) error_read_line(glob_fname_nuTable,ii+1);

    glob_map_suffixes[ii]=suf;
    glob_nu0_arr[ii]=nu0;
    glob_nuf_arr[ii]=nu1;
    glob_nu_arr[ii]=0.5*(nu0+nu1);
  }
  printf(" %d frequency channels found, labelled %03d to %03d\n\n",
	 glob_n_nu,glob_map_suffixes[0],glob_map_suffixes[glob_n_nu-1]);
}

static void read_rms(void)
{
  int ii;
  glob_rms_arr=(double *)my_malloc(glob_n_nu*sizeof(double));

#ifdef _RMS_FROM_MAPS
  glob_calculate_rms=1;
#else //_RMS_FROM_MAPS
  glob_calculate_rms=0;
#endif //_RMS_FROM_MAPS
  if(glob_calculate_rms) {
    for(ii=0;ii<glob_n_nu;ii++) {
      char fname[256];
      float *dum_map;
      long ns,jj;
      double mean=0,var=0;
      sprintf(fname,"%s_%03d.fits",glob_prefix_cosmo_in,glob_map_suffixes[ii]);
      dum_map=he_read_healpix_map(fname,&ns,0);
      if(ns!=glob_nside) report_error(1,"Shit\n");
      for(jj=0;jj<glob_nth;jj++) {
	double x=(double)(dum_map[glob_indices_in_mask[jj]]);
	mean+=x;
	var+=x*x;
      }
      mean/=glob_nth;
      glob_rms_arr[ii]=sqrt(var/glob_nth-mean*mean);
      free(dum_map);
    }

    FILE *fi=my_fopen(glob_fname_rms,"w");
    for(ii=0;ii<glob_n_nu;ii++) {
      fprintf(fi,"%lE %lE\n",glob_nu_arr[ii],glob_rms_arr[ii]);
    }
    fclose(fi);
  }
  else {
    FILE *fi=my_fopen(glob_fname_rms,"r");
    int nl=linecount(fi);
    rewind(fi);
    if(nl!=glob_n_nu) report_error(1,"Wrong number of frequency channels\n");
    for(ii=0;ii<glob_n_nu;ii++) {
      double nu,rms;
      int stat=fscanf(fi,"%lf %lf",&nu,&rms);
      if(stat!=2) error_read_line(glob_fname_rms,ii+1);
      glob_rms_arr[ii]=rms;
    }
    fclose(fi);
  }
  printf("\n");
}

void write_cleaned_maps()
{
  int ii;
  float *map_clean=(float *)my_calloc(12*glob_nside*glob_nside,sizeof(float));

  printf("Writing clean maps\n");
  for(ii=0;ii<glob_n_nu;ii++) {
    long jj;
    char fname[256];
    for(jj=0;jj<glob_nth;jj++)
      map_clean[glob_indices_in_mask[jj]]=(float)(glob_clean_maps[ii+jj*glob_n_nu]);
    sprintf(fname,"%s_clean_%03d.fits",glob_prefix_out,glob_map_suffixes[ii]);
    he_write_healpix_map(&map_clean,1,glob_nside,fname);
  }

  free(map_clean);
  printf("\n");
}

static void read_mask(void)
{
  long ii;

  glob_mask=he_read_healpix_map(glob_fname_mask,&glob_nside,0);
  glob_nth=0;
  
  printf(" * Resolution parameter Nside = %ld\n",glob_nside);

  for(ii=0;ii<12*glob_nside*glob_nside;ii++) {
    if((glob_mask[ii]<-0.01)||((glob_mask[ii]>0.01)&&(glob_mask[ii]<0.99))||(glob_mask[ii]>1.01))
      report_error(1,"Mask must be 0 or 1 : lE\n",glob_mask[ii]);
    if(glob_mask[ii]>0.99) {
      glob_mask[ii]=1;
      glob_nth++;
    }
    else glob_mask[ii]=0;
  }
  
  glob_fsky=(double)glob_nth/(12*glob_nside*glob_nside);
  printf(" * Sky fraction : %.2lf%%\n",glob_fsky*100);
  printf("\n");
  
  glob_indices_in_mask=(long *)my_malloc(glob_nth*sizeof(long));
  glob_nth=0;
  for(ii=0;ii<12*glob_nside*glob_nside;ii++) {
    if(glob_mask[ii]>0.99) {
      glob_indices_in_mask[glob_nth]=ii;
      glob_nth++;
    }
  }
}

static void read_params(char *fname_params)
{
  FILE *fi;
  int n_lin,ii;
  
  //Read parameters from file
  fi=my_fopen(fname_params,"r");
  n_lin=linecount(fi);
  rewind(fi);
  for(ii=0;ii<n_lin;ii++) {
    char s0[512],s1[64],s2[256];
    if(fgets(s0,sizeof(s0),fi)==NULL)
      error_read_line(fname_params,ii+1);
    if((s0[0]=='#')||(s0[0]=='\n')) continue;
    int sr=sscanf(s0,"%s %s",s1,s2);
    if(sr!=2)
      error_read_line(fname_params,ii+1);

    if(!strcmp(s1,"output_prefix="))
      sprintf(glob_prefix_out,"%s",s2);
    else if(!strcmp(s1,"input_prefix="))
      sprintf(glob_prefix_in,"%s",s2);
    else if(!strcmp(s1,"cosmo_input_prefix="))
      sprintf(glob_prefix_cosmo_in,"%s",s2);
    else if(!strcmp(s1,"fname_nutable="))
      sprintf(glob_fname_nuTable,"%s",s2);
    else if(!strcmp(s1,"fname_mask="))
      sprintf(glob_fname_mask,"%s",s2);
    else if(!strcmp(s1,"fname_rms="))
      sprintf(glob_fname_rms,"%s",s2);
    else if(!strcmp(s1,"method="))
      sprintf(glob_method,"%s",s2);
    else if(!strcmp(s1,"n_foregrounds="))
      glob_n_foregrounds=atoi(s2);
    else if(!strcmp(s1,"analyze_pk="))
      glob_analyze_pk=atoi(s2);
    else if(!strcmp(s1,"output_maps="))
      glob_output_maps=atoi(s2);
    else if(!strcmp(s1,"start_radial="))
      glob_start_radial=atoi(s2);
    else if(!strcmp(s1,"end_radial="))
      glob_end_radial=atoi(s2);    
    else if(!strcmp(s1,"n_bins_radial="))
      glob_n_bins_radial=atoi(s2);    
#ifdef _RMS_FROM_MAPS
    else if(!strcmp(s1,"calculate_rms="))
      glob_calculate_rms=atoi(s2);    
#endif //_RMS_FROM_MAPS
    else
      report_error(1,"Unknown parameter %s\n",s1);
  }
  fclose(fi);

  printf(" Run parameters:\n");
  printf("  * Output prefix: %s\n",glob_prefix_out);
  printf("  * Input prefix: %s\n",glob_prefix_in);
  printf("  * Input cosmo prefix: %s\n",glob_prefix_cosmo_in);
  printf("  * Frequency table file: %s\n",glob_fname_nuTable);
  printf("  * Mask file: %s\n",glob_fname_mask);
  printf("  * Inverse weights file %s\n",glob_fname_rms);
  if(!strcmp(glob_method,"pca"))
    printf("  * Method: PCA\n");
  else if(!strcmp(glob_method,"polog"))
    printf("  * Method: LOS fitting\n");
  else if(!strcmp(glob_method,"fastica"))
    printf("  * Method: ICA\n");
  else {
    printf("WTF??\n");
    exit(1);
  }
  printf("  * Number of FG dof: %d\n",glob_n_foregrounds);
  if(glob_analyze_pk)
    printf("  * %d independent bins for radial Pk\n",glob_n_bins_radial);
  if(glob_output_maps)
    printf("  * Clean maps written to %s_clean_XXX.fits\n",glob_prefix_out);
  printf("\n");
}

void read_data(char *fname_params)
{
  int ii;

  printf("Reading params\n");
  read_params(fname_params);
  printf("Reading frequency table\n");
  read_nu_arr();
  printf("Reading mask\n");
  read_mask();
  printf("Reading inverse weights\n");
  read_rms();

  printf("Reading maps\n");
  glob_data_maps=(double *)my_malloc(glob_nth*glob_n_nu*sizeof(double));
  glob_clean_maps=(double *)my_malloc(glob_nth*glob_n_nu*sizeof(double));
  for(ii=0;ii<glob_n_nu;ii++) {
    char fname[256];
    float *dum_map;
    long ns,jj;
    sprintf(fname,"%s_%03d.fits",glob_prefix_in,glob_map_suffixes[ii]);
#ifdef _VERBOSE
    printf(" Reading map %s\n",fname);
#endif //_VERBOSE
    dum_map=he_read_healpix_map(fname,&ns,0);
    if(ns!=glob_nside) report_error(1," Wrong resolution parameter : %d != %d\n",ns,glob_nside);
    for(jj=0;jj<glob_nth;jj++)
      glob_data_maps[ii+jj*glob_n_nu]=(double)(dum_map[glob_indices_in_mask[jj]]);
    free(dum_map);
  }
  printf("\n");
}

void postprocess_clean_maps(void)
{
  long ii;
  double *mean=(double *)my_calloc(glob_n_nu,sizeof(double));

  printf("Subtracting leftover mean from clean maps\n");

  for(ii=0;ii<glob_nth;ii++) {
    int jj;
    for(jj=0;jj<glob_n_nu;jj++)
      mean[jj]+=glob_clean_maps[jj+glob_n_nu*ii];
  }

  for(ii=0;ii<glob_n_nu;ii++)
    mean[ii]/=glob_nth;

  for(ii=0;ii<glob_nth;ii++) {
    int jj;
    for(jj=0;jj<glob_n_nu;jj++)
      glob_clean_maps[jj+glob_n_nu*ii]-=mean[jj];
  }

  free(mean);
  printf("\n");
}

void end_all(void)
{
  free(glob_data_maps);
  free(glob_clean_maps);
  free(glob_nu0_arr);
  free(glob_nuf_arr);
  free(glob_nu_arr);
  free(glob_map_suffixes);
  free(glob_indices_in_mask);
  free(glob_mask);
}
