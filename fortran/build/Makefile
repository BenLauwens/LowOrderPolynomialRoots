.SUFFIXES: .f90

EXE = fcsprp

F95 = ifort -O3 -xHost -fast -warn -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
F95 = gfortran -O3 -march=native -flto -Wall -fstack-arrays

VPATH = ../src

all : $(EXE)

.PHONY: clean
clean :
	rm *.mod *.o $(EXE)

run : $(EXE)
	sudo nice -1 ./$(EXE) $(CMD)

$(EXE) : main.o
	$(F95) -o $(EXE) *.o

main.o : fcsprp.o

%.o : %.f90
	$(F95) -c $<
