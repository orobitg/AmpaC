PROJECT = AmpaC

# Environment 
SRC = ./src
OBJ = ./obj
BIN = ./bin

GCC = gcc
MKDIR=mkdir
CP=cp
#GCCFLAGS := -I$(INC) -I$(MPIEXT)

# Add any options you like
GCCFLAGS += -g
GCCFLAGS += -O3
GCCFLAGS += -lm #-lpthread

OBJS =	main.o
SRCS = main.c
	  
vpath %.c $(SRC)
vpath %.o   $(OBJ)
vpath %     $(BIN)

all: build

build: $(PROJECT)

#dir:
#	mkdir -p $(BIN); mkdir -p $(OBJ);

clean:
	rm -f *.o

$(PROJECT): $(OBJS)
	$(GCC) $(GCCFLAGS) $(SRCS) -o $(PROJECT)

