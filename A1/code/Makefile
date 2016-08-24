####################################################################
#                                                                  #
#             Makefile for FIRE Solver Benchmarks                  #
#                                                                  #
####################################################################

CC = mpicc -icc
CFLAGS = -Wall -g -O1
LIBS = -lm $(PAPI_INC) $(PAPI_LIB)

LIBPOS=libpos.a
AR = ar
ARFLAGS = rv

SRCS = initialization.c compute_solution.c finalization.c util_read_files.c util_write_files.c
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: gccg 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

gccg: gccg.c $(LIBPOS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(LIBPOS) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

clean:
	rm -rf *.o gccg $(LIBPOS) binconv


binconv:
	$(CC) -o binconv text2bin.c
