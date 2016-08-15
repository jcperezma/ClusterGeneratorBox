FC90 = ifort -heap-arrays -traceback -fpp -DSHEAR_FLOW
FCFLAGS=  -i8 -I${MKLROOT}/include -mkl:sequential  
MKLROOT = /home/pec1065/intel/compilers_and_libraries_2016.1.150/linux/mkl
LIBS =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm  
TARGETS= clean ClusterGen
OBJSOC = EXCL_VOL.o
OBJMOD = $(wildcard *.mod)

ClusterGen: $(OBJSOC) main.f90 
		$(FC90) $(FCFLAGS) -o ClusterGenerator main.f90  $(LIBS) $(OBJSOC)
		
EXCL_VOL.o: ex_vol.f90
		$(FC90) -c ex_vol.f90 

		
clean: 
	rm -f *.o *.mod

