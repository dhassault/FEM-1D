# TP2 finite element method 1D
FC = gfortran
FFLAGS = -Wall
LDFLAGS = -Wall

SRCS = fem1d.f90
OBJS = $(SRCS:.f90=.o)
EXEC = $(SRCS:.f90=)
MOD = $(SRCS:.f90=.mod)
VPATH = libs/

all: $(EXEC)

%.o: %.f90
					$(FC) $(FFLAGS) -c $<


%: %.o
					$(FC) $(FFLAGS) -o $@ $^
					rm -f $(OBJS) *.o *.mod
fem1d: functions.o linear_system.o discretisation.o solver.o save_data.o estimation_error.o


clean :
					rm -f $(EXEC) *.o *.mod
