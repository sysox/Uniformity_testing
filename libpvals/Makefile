TARGET=libpvals
CPPFLAGS=
CFLAGS=-O3 -g -Wall -fPIC -D_GNU_SOURCE
LDLIBS=-lm -lgsl
CC=gcc

SOURCES=\
    dieharder/dieharder-kstest.c dieharder/dieharder.c \
    nist-sts/nist-sts-cephes.c nist-sts/nist-sts-assess.c \
    testu01/fbar.c testu01/fdist.c testu01/util.c \
    testu01/num.c testu01/num2.c testu01/fmass.c \
    testu01/tables.c testu01/mystr.c testu01/bitset.c \
    testu01/gofs.c testu01/gofw.c testu01/wdist.c testu01/testu01.c \
    utils.c

OBJECTS=$(SOURCES:.c=.o)

all: $(TARGET).a $(TARGET).so

$(TARGET).a: $(OBJECTS)
	ar rcs $@ $(OBJECTS)

$(TARGET).so: $(OBJECTS) libpvals.sym
	$(CC) $(CPPFLAGS) $(OBJECTS) -o $@ $(LDLIBS) -shared -Wl,--version-script=libpvals.sym

clean:
	rm -f *.o dieharder/*.o nist-sts/*.o testu01/*.o *~ core *.a *.so

.PHONY: clean
