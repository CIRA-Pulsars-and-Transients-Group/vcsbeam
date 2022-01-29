CC=g++
LDLIBS=-lvcsbeam -lcudart -lcufft -lmwalib
CFLAGS=-I/home/smcsweeney/include \
	   -I/opt/cuda/targets/x86_64-linux/include \
	   $(shell mpicc -showme:link)
LDFLAGS=-L/opt/cuda/targets/x86_64-linux/lib

TARGETS = applycal

all: $(TARGETS)

clean:
	$(RM) $(TARGETS)
