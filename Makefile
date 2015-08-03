#ViennaRNA package location
ViennaRNA=~/software/ViennaRNA-2.1.8/H
ViennaLIB=~/software/ViennaRNA-2.1.8/lib/libRNA.a

CPP = g++
CC = gcc
CPPFLAGS =  -O2 -g -Wall -std=c++0x -fexceptions -Wno-write-strings
CFLAGS =  -O2 -g -Wall -fexceptions -Wno-write-strings
LFLAGS = -O2 -fopenmp

OBJ = RNAlocmin_cmdline.o\
			barrier_tree.o\
			flood.o\
			hash_util.o\
			main.o\
			RNAlocmin.o\
			treeplot.o\
			globals.o\
			move_set_pk.o\
			pknots.o\
			findpath_pk.o\
			move_set_inside.o

DIRS = -I $(ViennaRNA)

LIBS = $(ViennaLIB)

all: $(OBJ)
	$(CPP) $(LFLAGS) $(DIRS) $(OBJ) $(LIBS) -o RNAlocmin
	rm -f $(OBJ)

RNAlocmin_cmdline.h RNAlocmin_cmdline.c: RNAlocmin.ggo
	gengetopt -i RNAlocmin.ggo

%.o: %.cpp
	$(CPP) $(CPPFLAGS) $(DIRS) -c $<

%.o: %.c
	$(CC) $(CFLAGS) $(DIRS) -c $<

clean:
	rm -f $(OBJ)
	rm -f RNAlocmin
	rm -f RNAlocmin_cmdline.c RNAlocmin_cmdline.h

