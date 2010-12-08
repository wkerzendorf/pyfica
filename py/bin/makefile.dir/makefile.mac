#### definitions

#FORTRAN = gfortran
FORTRAN = ifort
#FORTRAN = ~/intel/fc/10.1.015/bin/ifort
FFLAGS = -O3 -w -r8 -Vaxlib
#FFLAGS  = -w -ffixed-line-length-none
#FFLAGS = -e -w -vax -g -xtypemap=real:64
EXNAME = fica.exe

#### macros:

OBJS = fica.o ficop.o ficj.o fick.o ficl.o fici.o ficit.o fiym.o \
       ficg.o fiyn.o ficst.o ficp.o fics.o ficf.o fico.o fiyq.o ficr.o \
       fiyt.o fiyz.o consts.o sort1.o bnu1.o wair1.o wvac.o ficcl.o

#### targets:

all: $(OBJS)
	$(FORTRAN) $(FFLAGS) $(OBJS) -o $(EXNAME)

fica.o: fica.f
	$(FORTRAN) $(FFLAGS) -c fica.f

ficop.o: ficop.f
	$(FORTRAN) $(FFLAGS) -c ficop.f

ficj.o: ficj.f
	$(FORTRAN) $(FFLAGS) -c ficj.f

fick.o: fick.f
	$(FORTRAN) $(FFLAGS) -c fick.f

ficl.o: ficl.f
	$(FORTRAN) $(FFLAGS) -c ficl.f

fici.o: fici.f
	$(FORTRAN) $(FFLAGS) -c fici.f

ficit.o: ficit.f
	$(FORTRAN) $(FFLAGS) -c ficit.f

fiym.o: fiym.f
	$(FORTRAN) $(FFLAGS) -c fiym.f

ficg.o: ficg.f
	$(FORTRAN) $(FFLAGS) -c ficg.f

fiyn.o: fiyn.f
	$(FORTRAN) $(FFLAGS) -c fiyn.f

ficst.o: ficst.f
	$(FORTRAN) $(FFLAGS) -c ficst.f

ficp.o: ficp.f
	$(FORTRAN) $(FFLAGS) -c ficp.f

fics.o: fics.f
	$(FORTRAN) $(FFLAGS) -c fics.f

ficf.o: ficf.f
	$(FORTRAN) $(FFLAGS) -c ficf.f

fico.o: fico.f
	$(FORTRAN) $(FFLAGS) -c fico.f

fiyq.o: fiyq.f
	$(FORTRAN) $(FFLAGS) -c fiyq.f

ficr.o: ficr.f
	$(FORTRAN) $(FFLAGS) -c ficr.f

fiyt.o: fiyt.f
	$(FORTRAN) $(FFLAGS) -c fiyt.f

fiyz.o: fiyz.f
	$(FORTRAN) $(FFLAGS) -c fiyz.f

consts.o: consts.f
	$(FORTRAN) $(FFLAGS) -c consts.f

sort1.o: sort1.f
	$(FORTRAN) $(FFLAGS) -c sort1.f

bnu1.o: bnu1.f
	$(FORTRAN) $(FFLAGS) -c bnu1.f

wair1.o: wair1.f
	$(FORTRAN) $(FFLAGS) -c wair1.f

wvac.o: wvac.f
	$(FORTRAN) $(FFLAGS) -c wvac.f

ficcl.o: ficcl.f
	$(FORTRAN) $(FFLAGS) -c ficcl.f

#### shell commands

clean:;	rm -f $(OBJS) fiym.o fiyn.o core 

