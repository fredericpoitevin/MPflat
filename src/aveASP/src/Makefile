#
EXT=
#NAME=AquaSAXS
NAME=aveASP

DIR=../bin/
NAMEFUL=$(DIR)$(NAME)$(EXT)

FC = gfortran
FFLAGS = -c -u -C -O -fno-align-commons -ffixed-line-length-none 

LDFLAGS = -O

.f.o:
	$(FC) $(FFLAGS) $<

OBJECTS = \
$(NAME).o \
read_input.o cubature.o atomparam.o \
surftens.o \
svdfit.o svdcmp.o svbksb.o pythag.o \
surface.o vector.o

#OBJECTS = \
#$(NAME).o \
#read_input.o cubature.o atomparam.o \
#profilexp.o atomicsf.o solvterm.o \
#saxs.o svdfit.o svdcmp.o svbksb.o pythag.o \
#writeprof.o maputils.o surface.o vector.o

$(NAME) : $(OBJECTS)
	$(FC) -o $(NAMEFUL) $(LDFLAGS) $(OBJECTS)


all : $(OBJECTS)
	$(FC) -o $(NAMEFUL) $(LDFLAGS) $(OBJECTS)

clean :
	touch junk.o; rm -f *.o $(NAMEFUL)

$(OBJECTS) : param.h
