CC=clang -O0 -Wall -ggdb
LIBS=-lm

OBJS=vortices.o vec3_maths.o tangle.o vortex_utils.o

all: $(OBJS)
	$(CC) $(OBJS) -o vortices $(LIBS)

vorices.o: vortices.c
	$(CC) -c vortices.c

vortex_utils.o: vortex_utils.c tangle.o vec3_maths.o
	$(CC) -c vortex_utils.c

tangle.o: tangle.c vec3_maths.o vortex_constants.h
	$(CC) -c tangle.c

vec3_maths.o: vec3_maths.c
	$(CC) -c vec3_maths.c


clean:
	rm vortices vortices.o
