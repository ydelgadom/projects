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
#			make PAR=mu anal
#			make PAR=tau anal
#   clean: make clean
# 
###############################################################
SIZE:=8
PAR:=tau

CFLAGS := -O3 -Wall

ran : ranlxd.c
	gcc $(CFLAGS) -c ranlxd.c

gen: ranlxd.o gen.cpp
	g++ $(CFLAGS) -D$(PAR) -DSIZE=$(SIZE) -o ./bin/gen$(SIZE)_$(PAR).x gen.cpp ranlxd.o

anal: analysis.cpp
	g++ $(CFLAGS) -D$(PAR) -o ./bin/anal_$(PAR).x analysis.cpp

clean:
	rm -v *.o ./bin/*.x
