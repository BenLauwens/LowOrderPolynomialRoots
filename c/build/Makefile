.SUFFIXES: .f90

EXE = fcsprp

CC = icc -xHost -O3 -std=c11 -fast -I/usr/local/lib -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
CC = gcc -march=native -O3 -flto -Wall -std=c17
#CC = clang -march=native -O3 -flto -Wall -std=c17

INCLUDE = -lgsl -lgslcblas -lm

VPATH = ../src

DEPS = fcsprp.h

OBJ = main.o fcsprp.o

all : $(EXE)

.PHONY: clean
clean :
	rm *.o $(EXE)

run : $(EXE)
	sudo nice -1 ./$(EXE)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< 

$(EXE): $(OBJ)
	$(CC) $(INCLUDE) -o $@ $^ 
