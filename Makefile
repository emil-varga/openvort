#####################################################################
# Copyright (C) 2018 Emil Varga <varga.emil@gmail.com>
# 
# This file is part of OpenVort
# 
# OpenVort is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#####################################################################

CC=gcc -fopenmp
#CFLAGS=  -O2 -fsanitize=undefined -fsanitize=float-divide-by-zero -fsanitize=leak `pkg-config --cflags libconfig`
CFLAGS = -O2 `pkg-config --cflags libconfig`
WARNINGS = -Wall -Wextra -pedantic
DEBUG = -ggdb -D_DEBUG_
INCLUDE = include
BUILD = build

COMPILE = $(CC) $(CFLAGS) $(WARNINGS) $(DEBUG) -I$(INCLUDE)
LIBS=-lm `pkg-config --libs libconfig`

src = $(wildcard src/*.c)
obj = $(patsubst src/%.c, $(BUILD)/%.o, $(src))
dep = $(obj:.o=.d)

all: vortices

vortices: $(obj)
	$(COMPILE) $^ -o $@ $(LIBS)

-include $(dep)

$(BUILD)/%.d: src/%.c | $(BUILD)
	$(CPP) $(CFLAGS) $< -I$(INCLUDE) -MM -MT $(@:.d=.o) >$@

$(BUILD)/%.o: src/%.c | $(BUILD)
	$(COMPILE) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

.PHONY: clean
clean:
	rm -rf vortices $(BUILD) $(dep)
