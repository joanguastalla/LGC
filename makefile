# Prefix used to recognize rules in recipes
.RECIPEPREFIX= 

#Define C Compiler
CC=gcc -O3  

# Compiler flags:
# -g Add debuging information to the executable file
# -Wall turn on most of compiler warnings
CFLAGS=-Wall -pedantic

# Define Library paths to include
LFLAGS=-L/usr/lib/x86_64-linux-gnu 

# Define libraries to be used
LIBS= -lm -lblas -llapack 

#Link openmp to gcc compiler
OMPLINK=-fopenmp

#Object for compilation
OBJ=modelling_2D.o modelling_utils.o

# Name of executable file
MAIN=modelling

all: $(MAIN) 
	@echo Compiling executable $(MAIN)

%.o: %.c 
	$(CC) -c -o $@ $< $(CFLAGS)

$(MAIN): $(OBJ)  
	$(CC) -o $@ $^ $(LIBS) $(LFLAGS)



clean:
	rm -f *.o $(MAIN)
