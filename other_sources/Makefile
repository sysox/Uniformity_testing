TARGET=test-ks
CFLAGS=-O3 -g -Wall -fPIC
LDLIBS=-lm
CC=gcc

SOURCES=k.c main.c

OBJECTS=$(SOURCES:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean:
	rm -f *.o k.o *~ core *.a *.so $(TARGET)

.PHONY: clean
