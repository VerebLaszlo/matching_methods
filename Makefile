CC=gcc

include config.mk

BUILD_TYPE=debug# debug, normal, prod
ifeq (${BUILD_TYPE},debug)
	CFLAGS+=-Wall -Wextra -Wswitch-default -Wswitch-enum -g3 -DDEBUG
else
ifeq (${BUILD_TYPE},prod)
	CFLAGS+=-O3
else
	CFLAGS+=-Wall -Wextra -Wswitch-default -Wswitch-enum -g3 -ggdb3 -time
endif
endif

OBJ=util_math.o datatypes.o fitting.o
CFLAGS+=$(shell pkg-config --cflags lalinspiral)
LDFLAGS+=$(shell pkg-config --libs lalinspiral)

all: main wave

run: rmain rwave

main: main.c ${OBJ}
	${CC} -o main ${CFLAGS} main.c ${OBJ} ${LDFLAGS}
	@echo ''

simple-test: simple-test.c ${OBJ}
	${CC} -o simple-test ${CFLAGS} simple-test.c ${OBJ} ${LDFLAGS}

fitting.o: fitting.c fitting.h util_math.o datatypes.o
	${CC} -c ${CFLAGS} fitting.c ${LDFLAGS}
	@echo ''

%.o: %.c %.h
	${CC} -c ${CFLAGS} $<
	@echo ''

rmain: main
	./main `head -n 1 input.data` stat.txt
	@echo ''

wave: LALSQTPNWaveformTest.c
	${CC} -o wave ${CFLAGS} LALSQTPNWaveformTest.c ${LDFLAGS}
	@echo ''

rwave: wave
	./wave `head -n 1 input.data`
	@echo ''

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
	@echo 'clean     : deletes the object files'
	@echo 'cleanrun : deletes the exe files'
	@echo 'cleanall : invokes the "clean" and "cleanrun" commands'
	@echo 'help      : prints this message'
	@echo ''

.PHONY: all clean cleanall cleanrun help
