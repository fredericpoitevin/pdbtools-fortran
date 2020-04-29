#
EXT=.gfort
NAME=lrt

DIR=../bin/
NAMEFUL=$(DIR)$(NAME)$(EXT)

FC = gfortran
FFLAGS = -c -u -C -O
LDFLAGS = -O

.f.o:
	$(FC) $(FFLAGS) $<

OBJECTS = \
$(NAME).o \
svdcmp.o pythag.o

$(NAME) : $(OBJECTS)
	$(FC) -g -o $(NAMEFUL) $(LDFLAGS) $(OBJECTS)

all : $(OBJECTS)
	$(FC) -g -o $(NAMEFUL) $(LDFLAGS) $(OBJECTS)

clean :
	touch junk.o; rm -f *.o $(NAMEFUL)

$(OBJECTS) : lrt.prm
