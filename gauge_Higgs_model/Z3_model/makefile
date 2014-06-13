##############################################################
#
#		Makefile
#		
#		To compile:
#   random number generator (.o file): make ran
#		program to generate configurations:  
#			make NS=4 NT=4 PAR=BETA gen
#			make NS=4 NT=4 PAR=MU gen
#		program to compute observables: 
#			make PAR=MU anal
#			make PAR=BETA anal
#   clean: make clean
# 
###############################################################
NS:=4
NT:=4
PAR:=BETA

CFLAGS := -O3 -Wall

ran : ranlxd.c
	gcc $(CFLAGS) -c -o ranlxd.o ranlxd.c

gen: ranlxd.o gen.cpp
	g++ $(CFLAGS) -D$(PAR) -DNS=$(NS) -DNT=$(NT) -o ./bin/gen$(NS)x$(NT)_$(PAR).x gen.cpp ranlxd.o

anal: analysis.cpp
	g++ $(CFLAGS) -D$(PAR) -o ./bin/anal_$(PAR).x analysis.cpp 

clean:
	rm -v *.o ./bin/*.x
