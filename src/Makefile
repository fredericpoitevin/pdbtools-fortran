EXT=
NAME=PDBorganizer
NAME2=mindist
NAME3=fix_numb
NAME4=twist
NAME5=orient
NAME6=flip
NAME7=superimpose
NAME8=interp_movie
NAME9=pdb2ellips
#
DIR=../bin
NAMEFUL=$(DIR)/$(NAME)$(EXT)
NAMEFUL2=$(DIR)/$(NAME2)$(EXT)
NAMEFUL3=$(DIR)/$(NAME3)$(EXT)
NAMEFUL4=$(DIR)/$(NAME4)$(EXT)
NAMEFUL5=$(DIR)/$(NAME5)$(EXT)
NAMEFUL6=$(DIR)/$(NAME6)$(EXT)
NAMEFUL7=$(DIR)/$(NAME7)$(EXT)
NAMEFUL8=$(DIR)/$(NAME8)$(EXT)
NAMEFUL9=$(DIR)/$(NAME9)$(EXT)
#
FC = gfortran
FFLAGS = -c -O -fimplicit-none
LDFLAGS = -O

.f.o :
	$(FC) $(FFLAGS) $<
	
OBJECTS = \
$(NAME).o PDBorg_routines.o crd_handler.o

$(NAME) : $(OBJECTS)
	$(FC) -o $(NAMEFUL) $(LDFLAGS) $(OBJECTS)

OBJECTS2 = \
$(NAME2).o crd_handler.o

$(NAME2) : $(OBJECTS2)
	$(FC) -o $(NAMEFUL2) $(LDFLAGS) $(OBJECTS2)

OBJECTS3 = \
$(NAME3).o 

$(NAME3) : $(OBJECTS3)
	$(FC) -o $(NAMEFUL3) $(LDFLAGS) $(OBJECTS3)

OBJECTS4 = \
$(NAME4).o crd_handler.o

$(NAME4) : $(OBJECTS4)
	$(FC) -o $(NAMEFUL4) $(LDFLAGS) $(OBJECTS4)

OBJECTS5 = \
$(NAME5).o crd_handler.o jacobi.o

$(NAME5) : $(OBJECTS5)
	$(FC) -o $(NAMEFUL5) $(LDFLAGS) $(OBJECTS5)

OBJECTS6 = \
$(NAME6).o crd_handler.o jacobi.o

$(NAME6) : $(OBJECTS6)
	$(FC) -o $(NAMEFUL6) $(LDFLAGS) $(OBJECTS6)

OBJECTS7 = \
$(NAME7).o crd_handler.o compute.o jacobi.o

$(NAME7) : $(OBJECTS7)
	$(FC) -o $(NAMEFUL7) $(LDFLAGS) $(OBJECTS7)

OBJECTS8 = \
$(NAME8).o crd_handler.o

$(NAME8) : $(OBJECTS8)
	$(FC) -o $(NAMEFUL8) $(LDFLAGS) $(OBJECTS8)

OBJECTS9 = \
$(NAME9).o crd_handler.o

$(NAME9) : $(OBJECTS9)
	$(FC) -o $(NAMEFUL9) $(LDFLAGS) $(OBJECTS9)


clean:
	touch junk.o; rm -f *.o $(NAMEFUL) $(NAMEFUL2) $(NAMEFUL3) $(NAMEFUL4) $(NAMEFUL5) $(NAMEFUL6) $(NAMEFUL7) $(NAMEFUL8) $(NAMEFUL9)

$(OBJECTS) : param.h
$(OBJECTS2) : param.h
$(OBJECTS4) : param.h
$(OBJECTS5) : param.h
$(OBJECTS6) : param.h
$(OBJECTS7) : param.h
$(OBJECTS8) : param.h
$(OBJECTS9) : param.h
