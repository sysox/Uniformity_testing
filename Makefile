TARGET=Uniformity_testing
CPPFLAGS=
CFLAGS=-O3 -g -Wall
LDLIBS=-lm -lgsl libpvals/libpvals.a other_sources/libks.a
LDFLAGS=
CC=gcc

SOURCES=main.c

OBJECTS=$(SOURCES:.c=.o)

all: $(TARGET) other_sources/libks.a

libpvals/libpvals.a:
	make -C libpvals

other_sources/libks.a:
	make -C other_sources

$(TARGET): $(OBJECTS) libpvals/libpvals.a other_sources/libks.a
	$(CC) -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean:
	rm -f *.o *~ core $(TARGET)
	make -C libpvals clean
	make -C other_sources clean

.PHONY: clean
