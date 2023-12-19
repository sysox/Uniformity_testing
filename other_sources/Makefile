TARGET=other_sources
CPPFLAGS=
CFLAGS=-O3 -g -Wall -fPIC -D_GNU_SOURCE
LDLIBS=-lm -lgsl
CC=gcc

SOURCES=k.c

OBJECTS=$(SOURCES:.c=.o)

all: $(TARGET).a $(TARGET).so

$(TARGET).a: $(OBJECTS)
	ar rcs $@ $(OBJECTS)

$(TARGET).so: $(OBJECTS)
	$(CC) $(CPPFLAGS) $(OBJECTS) -o $@ $(LDLIBS) -shared -W

clean:
	rm -f *.o k.o *~ core *.a *.so

.PHONY: clean
