##############################################################
#
#		Makefile
#		
#		To compile:
#   random number generator (.o file): make ran
#		program to generate configurations:  
#			make SIZE=4 gen
#		program to compute observables: 
#			make anal
#   clean: make clean
# 
###############################################################
SIZE:=8

CFLAGS := -O3 -Wall
LIBS := -lgsl -lblas

ran : ranlxd.c
	gcc $(CFLAGS) -c -o ranlxd.o ranlxd.c

gen: ranlxd.o gen.cpp
	g++ $(CFLAGS) -DSIZE=$(SIZE) -o ./bin/gen$(SIZE).x gen.cpp ranlxd.o $(LIBS)

anal: analysis.cpp
	g++ $(CFLAGS) -o ./bin/anal.x analysis.cpp $(LIBS)

clean:
	rm -v *.o ./bin/*.x
