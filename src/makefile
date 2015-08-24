SHELL = /bin/sh

VPATH := ../src/
ODIR = odir/

MACHINE_NAME := $(shell hostname)   

GIT_SHA := $(shell git rev-parse HEAD)    
GIT_MOD := $(shell git diff-index --quiet HEAD; echo $$?)
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD )

GIT_MOD := $(strip $(subst 1,+,$(GIT_MOD)))
GIT_MOD := $(strip $(subst 0,,$(GIT_MOD)))

########################################################################
#                  Compiler Flags (Build Specific)                     #
########################################################################

ifeq ($(BUILDTYPE),DG_WASUPP)
  FC = ifort
  
  FFLAGS_VEC  =  -O3 -xHost -align array32byte -align rec32byte -vec-report6  -I$(ODIR) -c
  FFLAGS      =  -O3 -xHost -align array32byte -align rec32byte  -I$(ODIR) -c 
  FFLAGS_LINK =  -O3 -xHost -o 

  EXE_NAME = DG_WASUPP
  MESSAGE = Compling serial program with optimization flags
endif

ifeq ($(BUILDTYPE),DG_WASUPP_DEBUG)
  FC = ifort
  
  FFLAGS_VEC  =  -C -g -traceback -I$(ODIR) -c
  FFLAGS      =  -C -g -traceback -I$(ODIR) -c 
  FFLAGS_LINK =  -o 

  EXE_NAME = DG_WASUPP_debug
  MESSAGE = Compling serial program with debug flags
endif
########################################################################
#                  Library Links (Platform Specific)                   #
########################################################################

NAME := $(strip $(subst ., , $(MACHINE_NAME)))

ifeq ($(NAME),diaplo-admin)
  LAPACK_LIB = -llapack
endif

ifeq ($(filter crc nd edu,$(NAME)),crc nd edu)
  LAPACK_LIB = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a  $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
endif
########################################################################
#                               Objects                                #
########################################################################

objects = SIZES.o GLOBALS.o READ_DGINP.o PREP_DG.o DG_TIMESTEP.o MAIN.o \
          ERROR_CHECKS.o SOLUTION_IO.o version.o NONHYDRO.o

obj := $(patsubst %.o, $(ODIR)%.o,$(objects))

########################################################################
#                          Executable Targets                          #
########################################################################

.PHONY : $(ODIR) DG_WASUPP clean clobber

ifeq ($(MAKELEVEL),0)  # this occurs first to set the build specific info (compiler/flags)
                       # each target calls make again with appropriate settings

DG_WASUPP : 
	@echo "\n\n"
	$(MAKE) BUILDTYPE=DG_WASUPP DG_WASUPP
	
DG_WASUPP_DEBUG : 
	@echo "\n\n"
	$(MAKE) BUILDTYPE=DG_WASUPP_DEBUG DG_WASUPP
	
else                   # this is called the second time after the build specific stuff is set

#DG_WASUPP : $(ODIR) $(obj)
#	$(FC) $(FFLAGS_LINK) $(EXE_NAME) $(obj) $(METIS_LIB) $(LAPACK_LIB)
#	
#	mv $(EXE_NAME) ../work
	
DG_WASUPP : $(ODIR) version  $(obj)
	@echo  " \n\033[35m$(MESSAGE) on $(FULLNAME) \033[0m \n"
	@echo "Git SHA: $(GIT_SHA)"
	@echo "Modified: $(GIT_MOD) \n"
	$(FC) $(FFLAGS_LINK) $(EXE_NAME) $(obj) $(METIS_LIB) $(LAPACK_LIB)
	
	mv $(EXE_NAME) ../work	
	
endif

########################################################################
#                            Set-up Targets                            #
########################################################################

$(ODIR) :
	mkdir -p $@

version:
	@sed -i '9 c \      PRINT*, "  Branch: $(GIT_BRANCH)" ' $(VPATH)version.F90	
	@sed -i '10 c \      PRINT*, "  SHA: $(GIT_SHA) $(GIT_MOD)" ' $(VPATH)version.F90
	@sed -i '16 c \      WRITE(FORT16,*) "  Branch: $(GIT_BRANCH)" ' $(VPATH)version.F90	
	@sed -i '17 c \      WRITE(FORT16,*) "  SHA: $(GIT_SHA) $(GIT_MOD)" ' $(VPATH)version.F90

########################################################################
#                            Clean Targets                             #
########################################################################	
	
clean : 
	@echo "\n"
	rm -r odir
	
clobber :
	@echo "\n\n"
	rm -r odir
#	rm DG_WASUPP
	
########################################################################
#                            File Targets                              #
########################################################################

# $< = first prerequisite/depenency
# $@ = name of the target

$(ODIR)GLOBALS.o : GLOBALS.F90
	$(FC) $(FFLAGS) $< -o $@
	mv globals.mod $(ODIR)	

$(ODIR)SIZES.o : SIZES.F90
	$(FC) $(FFLAGS) $< -o $@
	mv sizes.mod $(ODIR)
	
$(ODIR)READ_DGINP.o : READ_DGINP.F90 SIZES.F90 GLOBALS.F90
	$(FC) $(FFLAGS) $< -o $@
	mv read_dginp.mod $(ODIR)
	
$(ODIR)version.o : version.F90
	$(FC) $(FFLAGS) $< -o $@	

$(ODIR)MAIN.o : MAIN.F90 GLOBALS.F90 SIZES.F90
	$(FC) $(FFLAGS) $< -o $@
	
$(ODIR)PREP_DG.o : PREP_DG.F90 GLOBALS.F90 SIZES.F90 READ_DGINP.F90
	$(FC) $(FFLAGS) $< -o $@
	
$(ODIR)NONHYDRO.o : NONHYDRO.F90 GLOBALS.F90 SIZES.F90 READ_DGINP.F90
	$(FC) $(FFLAGS) $< -o $@
	
$(ODIR)DG_TIMESTEP.o : DG_TIMESTEP.F90 GLOBALS.F90 SIZES.F90 READ_DGINP.F90 NONHYDRO.F90
	$(FC) $(FFLAGS) $< -o $@
	
$(ODIR)ERROR_CHECKS.o : ERROR_CHECKS.F90 GLOBALS.F90 SIZES.F90 READ_DGINP.F90
	$(FC) $(FFLAGS) $< -o $@
	
$(ODIR)SOLUTION_IO.o : SOLUTION_IO.F90 GLOBALS.F90 SIZES.F90 READ_DGINP.F90
	$(FC) $(FFLAGS) $< -o $@





