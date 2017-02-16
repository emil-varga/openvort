CC=clang -O0 -Wall -ggdb
LIBS=-lm

OBJS=vortices.o vec3_maths.o tangle.o vortex_utils.o util.o

all: $(OBJS)
	$(CC) $(OBJS) -o vortices $(LIBS)

vorices.o: vortices.c vec3_math.h tangle.h vortex_utils.h
	$(CC) -c vortices.c

vortex_utils.o: vortex_utils.c vortex_utils.h
	$(CC) -c vortex_utils.c

tangle.o: tangle.c tangle.h vortex_constants.h
	$(CC) -c tangle.c

vec3_maths.o: vec3_maths.c vec3_maths.h
	$(CC) -c vec3_maths.c

util.o: util.c
	$(CC) -c util.c

clean:
	rm vortices vortices.o
