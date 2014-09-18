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

  printf("|--------------|\n");
  printf("|     FG_RM    |\n");
  printf("|--------------|\n\n");

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

  printf("|--------------|\n");
  printf("|     Done!    |\n");
  printf("|--------------|\n");

  return 0;
}
