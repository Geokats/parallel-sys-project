# the compiler
CC = mpicc

# compiler flags
CFLAGS  = -O3

# the build target executable:
TARGET = mpi_heat2Dn

all: $(TARGET).x

$(TARGET).x: $(TARGET).c
	module load openmpi mpiP
	$(CC) $(CFLAGS) $(TARGET).c -o $(TARGET).x

clean:
	$(RM) $(TARGET).x
