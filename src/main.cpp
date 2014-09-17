#include "common.h"

void fgrm(void)
{
  if(!strcmp(glob_method,"pca"))
    do_pca(glob_n_foregrounds);
  else if(!strcmp(glob_method,"polog"))
    do_polog(glob_n_foregrounds);
  else if(!strcmp(glob_method,"fastica"))
    do_fastica(glob_n_foregrounds);
}

int main(int argc,char **argv)
{
  char fname_params[256];
  if(argc!=2) {
    fprintf(stderr,"Usage: fgrm param_file\n");
    exit(1);
  }
  sprintf(fname_params,"%s",argv[1]);

  setbuf(stdout,NULL);
  read_data(fname_params);
  fgrm();
  if(glob_analyze_pk) {
    do_analysis(glob_start_radial,glob_end_radial,
		glob_n_bins_radial);
  }
  if(glob_output_maps) {
    write_cleaned_maps();
  }
  end_all();

  return 0;
}
