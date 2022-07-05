FORTRAN=mpif90
OPTS= -O2 -Wall -Wtabs --free-form
LD=
LDOPTS=
EXENAME=nematic-membrane
OBJS=module_datastruct.o  module_rerun.o module_dataformat.o module_initialize_system.o module_linklist_calc.o  module_curvcalc.o module_energy_calculations.o module_mcsmoves.o  module_compute_analytic_measures.o maincode.o
$(EXENAME):$(OBJS)
	 $(FORTRAN) -o $(EXENAME) $(LD) $(LDOPTS) $(OPTS) $(OBJS)
.f.o:
	$(FORTRAN) $(OPTS) -c $<
clean:
	rm -fv $(EXENAME) $(OBJS) *.mod *.dat *.vtu *.jvx startdet-* rundir-*/*.vtu rundir-*/*.jvx rundir*/startdet-*.dat rundir*/*.jvx~ rundir*/*.dat
