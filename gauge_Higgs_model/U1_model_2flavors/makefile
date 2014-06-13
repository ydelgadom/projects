##############################################################
#
#		Makefile
#		
#		To compile:
#   random number generator (.o file): make ran
#		program to generate configurations:  
#			make NS=4 NT=8 gen
#		program to compute observables: 
#			make anal
#   clean: make clean
# 
###############################################################
NS:=4
NT:=4

CFLAGS:=-O3 -Wall
LIBS:=-lgsl -lblas
#-lboost
INCLUDE:=-I/hosts/nashome/ydd/lattice/libs/ -L/hosts/nashome/ydd/lattice/libs/

## compile rule for random generator
ranlxd.o : ranlxd.c
	gcc $(CFLAGS) -c ranlxd.c

gen: ranlxd.o gen.cpp
	mpicxx $(CFLAGS) -DNS=$(NS) -DNT=$(NT) -o bin/gen$(NS)x$(NT).x gen.cpp ranlxd.o $(LIBS) $(INCLUDE)

anal: analysis.cpp
	g++ $(CFLAGS) -o bin/anal.x analysis.cpp $(LIBS) $(INCLUDE)

clean:
	rm -v *.o bin/*.x
