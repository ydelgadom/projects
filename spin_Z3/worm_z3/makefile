##############################################################
#
#		Makefile
#		
#		To compile:
#   random number generator (.o file): make ran
#		program to generate configurations with CLOSED WORMS:  
#			make SIZE=4 PAR=TAU gen
#			make SIZE=4 PAR=MU gen
#			make SIZE=4 PAR=KAPPA gen
#		program to generate configurations with OPEN WORMS:  
#			make SIZE=4 PAR=TAU genopen
#		program to compute observables: 
#			make PAR=mu anal
#			make PAR=tau anal
#			make PAR=kappa anal
#   clean: make clean
# 
###############################################################
SIZE:=8
PAR:=TAU

CFLAGS := -O3 -Wall

ran : ranlxd.c
	gcc $(CFLAGS) -c -o ranlxd.o ranlxd.c

gen: ranlxd.o gen.cpp
	g++ $(CFLAGS) -DSIZE=$(SIZE) -D$(PAR) -o ./bin/gen$(SIZE)_$(PAR).x gen.cpp ranlxd.o

genopen: ranlxd.o gen.cpp
	g++ $(CFLAGS) -DSIZE=$(SIZE) -D$(PAR) -DOPEN -o ./bin/gen$(SIZE)_open_$(PAR).x gen.cpp ranlxd.o

anal: analysis.cpp
	g++ $(CFLAGS) -D$(PAR) -o ./bin/anal_$(PAR).x analysis.cpp

clean:
	rm -v *.o ./bin/*.x
