.phony:  all PROGRAM clean

# Choose the compiler and set compiler options
OMP_NUM_THREADS      := 1024#$(shell nproc)
export OMP_NUM_THREADS

#OMP_STACKSIZE        := 1024M
#export OMP_STACKSIZE

FCOMP     =  gfortran
FOPTS     =  -cpp -fopenmp -O3  -fexternal-blas  -march=native -w -fallow-argument-mismatch
# FOPTS    +=  -freal-8-real-10 
#FOPTS    +=  -fdefault-real-8

# specify BLAS and LAPACK library
LDOPT     =  -lblas -llapack -flto
#LDOPT     = OpenBLAS/libopenblas.a
#LDOPT     =  FLAME/libflame.a BLIS/libblis.a

FCOMPQUAD = $(FCOMP)
FC        = $(FCOMP)
FFLAGS    = $(FOPTS)

CC        = gcc
COPTS     = -O3  -I./include

export FC FFLAGS FCOMP FOPTS


# Set the list of programs to compile

PROGRAMS = experiment1 experiment2 experiment3 experiment4 eperiment5 experiment6 experiment_boxes test_levin2d_lege test_levin2d_final test_levin test_levin2d test_adapgauss2d


# Compile all of the test programs and the library
all	      	             : clean $(PROGRAMS) 


# List the dependencies for each module's test program


EXPERIMENT_BOXES_FILES         = utils.o                                                  \
                                 linalg0.o                                                \
                                 chebyshev.o                                              \
                                 chebpw.o                                                 \
                                 legendre.o                                               \
                                 levin.o                                                  \
                                 appell_nop.o                                             \
                                 riccati_nop.o                                            \
                                 besphase.o                                               \
                                 adapgauss2d.o                                            \
                                 levin2d.o

EXPERIMENT_BOX_FILES         = utils.o                                                    \
                                 linalg0.o                                                \
                                 chebyshev.o                                              \
                                 chebpw.o                                                 \
                                 legendre.o                                               \
                                 levin.o                                                  \
                                 appell_nop.o                                             \
                                 riccati_nop.o                                            \
                                 besphase.o                                               \
                                 adapgauss2d.o                                            \
                                 levin2d.o                                                
                                 
EXPERIMENT6_FILES              = utils.o                                                  \
                                 linalg0.o                                                \
                                 chebyshev.o                                              \
                                 chebpw.o                                                 \
                                 legendre.o                                               \
                                 levin.o                                                  \
                                 adapgauss2d.o                                            \
                                 levin2d.o

EXPERIMENT5_FILES              = utils.o                                                  \
                                 linalg0.o                                                \
                                 chebyshev.o                                              \
                                 chebpw.o                                                 \
                                 legendre.o                                               \
                                 levin.o                                                  \
                                 adapgauss2d.o                                            \
                                 levin2d.o

EXPERIMENT4_FILES              = utils.o                                                  \
                                 linalg0.o                                                \
                                 chebyshev.o                                              \
                                 chebpw.o                                                 \
                                 legendre.o                                               \
                                 levin.o                                                  \
                                 appell_nop.o                                             \
                                 riccati_nop.o                                            \
                                 besphase.o                                               \
                                 adapgauss2d.o                                            \
                                 levin2d.o                                          \

EXPERIMENT3_FILES              = utils.o                                                  \
                                 linalg0.o                                                \
                                 chebyshev.o                                              \
                                 chebpw.o                                                 \
                                 legendre.o                                               \
                                 levin.o                                                  \
                                 adapgauss2d.o                                            \
                                 levin2d.o

EXPERIMENT2_FILES              = utils.o                                                  \
                                 linalg0.o                                                \
                                 chebyshev.o                                              \
                                 chebpw.o                                                 \
                                 errfun.o                                                 \
                                 levin.o                                                  \
                                 chebpw.o                                                 \
                                 levin2d.o

EXPERIMENT1_FILES              = utils.o                                                  \
                                 linalg0.o                                                \
                                 chebyshev.o                                              \
                                 chebpw.o                                                 \
                                 errfun.o                                                 \
                                 levin.o                                                  \
                                 levin2d.o

ADAPGAUSS2D_FILES              = utils.o                                                  \
                                 legendre.o                                               \
                                 adapgauss2d.o                                            

LEVIN2D_LEGE_FILES             = utils.o                                                  \
                                 linalg0.o                                                \
                                 chebyshev.o                                              \
                                 legendre.o                                               \
                                 chebpw.o                                                 \
                                 levin.o                                                  \
                                 levin_lege.o                                             \
                                 levin2d_final.o                                          \
                                 levin2d_lege.o

LEVIN2D_FILES                  = utils.o                                                  \
                                 linalg0.o                                                \
                                 chebyshev.o                                              \
                                 legendre.o                                               \
                                 chebpw.o                                                 \
                                 levin.o                                                  \
                                 levin2d.o                                                \
                                 adapgauss2d.o 


LEVIN2D_FINAL_FILES            = utils.o                                                  \
                                 linalg0.o                                                \
                                 chebyshev.o                                              \
                                 chebpw.o                                                 \
                                 levin.o                                                  \
                                 levin2d_final.o



LEVIN_FILES                    = utils.o                                                  \
                                 linalg0.o                                                \
                                 chebyshev.o                                              \
                                 legendre.o                                               \
                                 chebpw.o                                                 \
                                 adapgauss.o                                              \
                                 levin.o

AMVW/amvw.a                    : 
	$(MAKE) -C AMVW

experiment_boxes.o            : $(EXPERIMENT_BOX_FILES) experiment_boxes.f90
experiment_boxes              : $(EXPERIMENT_BOX_FILES) experiment_boxes.o

experiment6.o                  : $(EXPERIMENT6_FILES) experiment6.f90
experiment6                    : $(EXPERIMENT6_FILES) experiment6.o

experiment5.o                  : $(EXPERIMENT5_FILES) experiment5.f90
experiment5                    : $(EXPERIMENT5_FILES) experiment5.o

experiment4.o                  : $(EXPERIMENT4_FILES) experiment4.f90
experiment4                    : $(EXPERIMENT4_FILES) experiment4.o

experiment3.o                  : $(EXPERIMENT3_FILES) experiment3.f90
experiment3                    : $(EXPERIMENT3_FILES) experiment3.o

experiment2.o                  : $(EXPERIMENT2_FILES) experiment2.f90
experiment2                    : $(EXPERIMENT2_FILES) experiment2.o

experiment1.o                  : $(EXPERIMENT1_FILES) experiment1.f90
experiment1                    : $(EXPERIMENT1_FILES) experiment1.o

test_adapgauss2d.o             : $(ADAPGAUSS2D_FILES) test_adapgauss2d.f90
test_adapgauss2d               : $(ADAPGAUSS2D_FILES) test_adapgauss2d.o

test_levin2d_lege.o           : $(LEVIN2D_LEGE_FILES) test_levin2d_lege.f90
test_levin2d_lege             : $(LEVIN2D_LEGE_FILES) test_levin2d_lege.o

test_levin2d_final.o           : $(LEVIN2D_FINAL_FILES) test_levin2d_final.f90
test_levin2d_final             : $(LEVIN2D_FINAL_FILES) test_levin2d_final.o

test_levin2d.o                 : $(LEVIN2D_FILES) test_levin2d.f90
test_levin2d                   : $(LEVIN2D_FILES) test_levin2d.o

test_levin.o                   : $(LEVIN_FILES) test_levin.f90
test_levin                     : $(LEVIN_FILES) test_levin.o


# Setup the general compilation rules

%		: %.o
	$(FCOMP) $(FOPTS) -o $@ $^ $(LDOPT)
	@echo  
	@echo 
	@echo "---------[ $@     ]--------------------------------------------"
	@echo 
	@./$@
	@echo 
	@echo "--------------------------------------------------------------------------"
	@echo 


%.o		: %.f90
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.f
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.cpp
	$(CPPCOMP) -c $(CPPOPTS)  $<

%.o		: %.c
	$(CC) -c $(COPTS)  $<

.PHONY: clean

clean:
	rm -f *.o *.mod *~ fort.* a.out *.a *.png
	rm -f $(PROGRAMS)
	rm -f *.dat


