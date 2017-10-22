CC=gcc -O0 -Wall -ggdb
CFLAGS= `gsl-config --cflags` -fsanitize=undefined -fsanitize=float-divide-by-zero -fsanitize=leak
LIBS=-lm `gsl-config --libs`

SRCS=tangle.c  util.c  vec3_maths.c  vortex_utils.c  vortices.c vortex_dynamics.c
OBJS=vortices.o vec3_maths.o tangle.o vortex_utils.o util.o vortex_dynamics.o

all: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o vortices $(LIBS)

depend: .depend

.depend: $(SRCS)
	echo $(SRCS)
	rm -f ./.depend
	$(CC) $(CFLAGS) -MM $^ > ./.depend;

include .depend

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm vortices *.o
