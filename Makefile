########## User-definable stuff ##########
#
###Compiler and compilation options
COMP_C = gcc
COMP_CPP = g++
OPTIONS = -Wall -O3 -fopenmp
#
### Behavioural flags
#Output more information (only necessary for debugging)
DEFINEFLAGS += -D_VERBOSE
#
###Path to libraries and headers
###If two or more of the dependencies reside in the same paths, only
###one instance is necessary.
#GSL
GSL_INC = -I/home/damonge/include
GSL_LIB = -L/home/damonge/lib
#FFTW
FFTW_INC = 
FFTW_LIB = 
#HEALPix
HPIX_INC = 
HPIX_LIB = 
#LIBSHARP
LSHT_INC =
LSHT_LIB =
#cfitsio
FITS_INC = 
FITS_LIB = 
#IT++
ITPP_INC = 
ITPP_LIB = 
#
########## End of user-definable ##########

ifeq ($(strip $(USE_SINGLE_PRECISION)),yes)
DEFINEFLAGS += -D_SPREC
LIB_FFTW = -lfftw3f_omp -lfftw3f
else
LIB_FFTW = -lfftw3_omp -lfftw3
endif

OPTIONS += $(DEFINEFLAGS) -D_SPREC

INC_ALL = -I./src $(GSL_INC) $(FFTW_INC) $(HPIX_INC) $(FITS_INC) $(LSHT_INC) $(ITPP_INC)
LIB_DIRS = $(GSL_LIB) $(FFTW_LIB) $(HPIX_LIB) $(FITS_LIB) $(LSHT_LIB) $(ITPP_LIB)
LIB_ALL = $(LIB_DIRS) -lgsl -lgslcblas -lfftw3 -lchealpix -lsharp -lfftpack -lc_utils -lcfitsio -litpp -lm

COMMONO = src/common.o
EHPIXO = src/healpix_extra.o
PKO = src/pk_analysis.o
PCAO = src/pca.o
POLOGO = src/polog.o
FASTICO = src/fastica.o
OFILES= $(COMMONO) $(EHPIXO) $(PKO) $(PCAO) $(POLOGO) $(FASTICO)

MAIN = src/main.cpp

EXE = fgrm

default : $(EXE)

%.o : %.c
	$(COMP_C) $(OPTIONS) $(INC_ALL) -c $< -o $@

%.o : %.cpp
	$(COMP_CPP) -Wno-write-strings $(OPTIONS) $(INC_ALL) -c $< -o $@

$(EXE) : $(OFILES)
	$(COMP_CPP) -Wno-write-strings $(OPTIONS) $(INC_ALL) $(OFILES) $(MAIN) -o $(EXE) $(LIB_ALL)

clean :
	rm -f src/*.o

cleaner : 
	rm -f *~ src/*.o src/*~ $(EXE)
