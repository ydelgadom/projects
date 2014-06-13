##############################################################
#
#		Makefile
#		
#		To compile:
#   random number generator (.o file): make ran
#		program to generate the configs.:  make SIZE=4 gen
#		program to do the analysis: make anal
#   clean: make clean
# 
#		Example:
#		make SIZE=16 gen
#		make anal	
#
###############################################################
SIZE:=8

CFLAGS := -O3 -Wall

ran : ranlxd.c
	gcc $(CFLAGS) -c ranlxd.c

gen: ranlxd.o local_metro_SU3.cpp
	g++ $(CFLAGS) -DSIZE=$(SIZE) -o ./bin/gen$(SIZE).x local_metro_SU3.cpp ranlxd.o

anal: local_metro_analyze.cpp
	g++ $(CFLAGS) -o ./bin/anal.x local_metro_analyze.cpp

clean:
	rm -v *.o ./bin/*.x
