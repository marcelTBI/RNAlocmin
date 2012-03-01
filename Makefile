#ViennaRNA package location
ViennaRNA=~/software/include/ViennaRNA
ViennaLIB=~/software/lib/libRNA.a

CPP = g++
CC = gcc
CFLAGS =  -O2 -g -Wall -fexceptions -Wno-write-strings
LFLAGS = -O2 -fopenmp

OBJ = RNAlocmin_cmdline.o\
			barrier_tree.o\
			flood.o\
			hash_util.o\
			main.o\
			move_set.o\
			treeplot.o\

DIRS =  -I $(ViennaRNA)\
				-I /usr/lib/gcc/x86_64-redhat-linux/4.6.1/include

LIBS = $(ViennaLIB)

all: $(OBJ)
	$(CPP) $(LFLAGS) $(DIRS) $(OBJ) $(LIBS) -o RNAlocmin
	rm -f $(OBJ)

RNAlocmin_cmdline.h RNAlocmin_cmdline.c: RNAlocmin.ggo
	gengetopt -i RNAlocmin.ggo

%.o: %.cpp
	$(CPP) $(CFLAGS) $(DIRS) -c $<

%.o: %.c
	$(CC) $(CFLAGS) $(DIRS) -c $<

clean:
	rm -f $(OBJ)
	rm -f RNAlocmin
	rm -f RNAlocmin_cmdline.c RNAlocmin_cmdline.h

