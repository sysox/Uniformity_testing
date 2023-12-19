TARGET=Uniformity_testing
CPPFLAGS=
CFLAGS=-O3 -g
LDLIBS=-lm -lgsl other_sources/KS.a
LDFLAGS=
CC=gcc

SOURCES=main.c

OBJECTS=$(SOURCES:.c=.o)

all: $(TARGET)

other_sources/KS.a:
	make -C libKS

$(TARGET): $(OBJECTS) other_sources/KS.a
	$(CC) -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean:
	rm -f *.o *~ core $(TARGET)
	make -C libKS clean

.PHONY: clean
