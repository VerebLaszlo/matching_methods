#!/usr/bin/colormake
############################################################################

CC=gcc
include config.mk

BUILD_TYPE=debug# debug, normal, prod
ifeq (${BUILD_TYPE},debug)
	CFLAGS+=-Wall -Wextra -Wswitch-default -Wswitch-enum -g3
	DEBUG=1
else
ifeq (${BUILD_TYPE},prod)
	CFLAGS+=-O3
	DEBUG=0
else
	CFLAGS+=-Wall -Wextra -Wswitch-default -Wswitch-enum -g3 -ggdb3 -time
	DEBUG=0
endif
endif

OBJ=util_math.o datatypes.o fitting.o
CFLAGS+=`pkg-config --cflags lalinspiral`
LIB=-lm `pkg-config --libs lalinspiral`

all: main wave

main: main.c ${OBJ}
	${CC} -o main ${CFLAGS} main.c ${OBJ} ${LIB}

fitting.o: fitting.c fitting.h util_math.o datatypes.o
	${CC} -c ${CFLAGS} fitting.c ${LIB} util_math.o datatypes.o

datatypes.o: datatypes.c datatypes.h
	${CC} -c ${CFLAGS} datatypes.c

util_math.o: util_math.c util_math.h
	${CC} -c ${CFLAGS} util_math.c -lm

run: wave
	./wave `head -n 1 input.data` own.out `tail -n 1 input.data`

wave: LALSQTPNWaveformTest.c
	${CC} -o SQT ${CFLAGS} LALSQTPNWaveformTest.c ${LIB}

clean:
	rm -rf *.o *.out *.b
	@echo ''

cleanrun:
	rm -rf main wave
	@echo ''

cleanall:
	make clean
	make cleanrun
	@echo ''

help :
	@echo 'all       : makes everything'
	@echo 'lal       : makes just the LALSTPNWaveform.c part'
	@echo 'own       : makes the whole own code'
	@echo 'clean     : deletes the object files'
	@echo 'cleanrun : deletes the exe files'
	@echo 'cleanall : invokes the "clean" and "cleanrun" commands'
	@echo 'run       : runs the two programs'
	@echo 'help      : prints this message'
	@echo ''

.PHONY: all clean cleanall cleanrun run help
