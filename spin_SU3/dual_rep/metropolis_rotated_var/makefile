##############################################################
#
#		Makefile
#		
#		To compile:
#   random number generator (.o file): make ran
#		program to generate configurations:  
#			make SIZE=4 PAR=TAU gen
#			make SIZE=4 PAR=MU gen
#			make SIZE=4 PAR=KAPPA0 gen
#		program to compute observables: 
#			make PAR=MU anal
#			make PAR=TAU anal
#			make PAR=KAPPA0 anal
#   clean: make clean
# 
###############################################################
SIZE:=8
PAR:=tau

CFLAGS := -O3 -Wall

ran : ranlxd.c
	gcc $(CFLAGS) -c ranlxd.c

gen: ranlxd.o gen.cpp
	g++ $(CFLAGS) -D$(PAR) -DSIZE=$(SIZE) -o ./bin/gen$(SIZE)_rot_$(PAR).x gen.cpp ranlxd.o

anal: analysis.cpp
	g++ $(CFLAGS) -D$(PAR) -o ./bin/anal_rot_$(PAR).x analysis.cpp

clean:
	rm -v *.o ./bin/*.x
