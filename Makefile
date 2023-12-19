TARGET=Uniformity_testing
CPPFLAGS=
CFLAGS=-O3 -g -Wall
LDLIBS=-lm -lgsl libpvals/libpvals.a
LDFLAGS=
CC=gcc

SOURCES=main.c

OBJECTS=$(SOURCES:.c=.o)

all: $(TARGET)

libpvals/libpvals.a:
	make -C libpvals
$(TARGET): $(OBJECTS) libpvals/libpvals.a
	$(CC) -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean:
	rm -f *.o *~ core $(TARGET)
	make -C libpvals clean

.PHONY: clean
