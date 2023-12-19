LIB=libks
TARGET=test-ks
CFLAGS=-O3 -g -Wall -fPIC
LDLIBS=-lm $(LIB).a
CC=gcc

TARGET_SOURCES=main.c
SOURCES=k.c

TARGET_OBJECTS=$(TARGET_SOURCES:.c=.o)
OBJECTS=$(SOURCES:.c=.o)

all: $(LIB).a $(LIB).so $(TARGET)

$(LIB).a: $(OBJECTS)
	ar rcs $@ $(OBJECTS)

$(LIB).so: $(OBJECTS)
	$(CC) $(CPPFLAGS) $(OBJECTS) -o $@ $(LDLIBS) -shared -W

$(TARGET): $(TARGET_OBJECTS) $(LIB).a
	$(CC) -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean:
	rm -f *.o k.o *~ core *.a *.so $(TARGET)

.PHONY: clean
