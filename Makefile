CC=gcc
CFLAGS=  -O0 `gsl-config --cflags` -fsanitize=undefined -fsanitize=float-divide-by-zero -fsanitize=leak
DEBUG = -Wall -ggdb -D_DEBUG_
INCLUDE = include

COMPILE = $(CC) $(CFLAGS) $(DEBUG) -I$(INCLUDE)
LIBS=-lm `gsl-config --libs`

src = $(wildcard src/*.c)
obj = $(src:.c=.o)
dep = $(obj:.o=.d)

all: vortices

vortices: $(obj)
	$(COMPILE) $^ -o $@ $(LIBS)

include $(dep)

%.d: %.c
	$(CPP) $(CFLAGS) $< -I$(INCLUDE) -MM -MT $(@:.d=.o) >$@

%.o: %.c
	$(COMPILE) -c $< -o $@

.PHONY: clean
clean:
	rm vortices $(obj) $(dep)
