PREFIX = /home/ryand/bin #This can be changed
INCLUDE_DIRS = htslib
CC = gcc
OPTS = -Wall -g -O3

OBJS = alignmentHeap.o countMinSketch.o graph.o murmur3.o TargetCreator.o realigner.o SemiGlobal.o threads.o TargetCreator_main.o realigner_main.o
VERSION = 0.1.0

#If we're building from a git repo, then append the most recent tag
ifneq "$(wildcard .git)" ""
VERSION := $(VERSION)-$(shell git describe --always --dirty)
endif

.PHONY: all clean htslib install clean-all version.h

.SUFFIXES:.c .o

all: MethIndelRealigner

.c.o:
	$(CC) -c $(OPTS) -I$(INCLUDE_DIRS) $< -o $@

version.h:
	echo '#define VERSION "$(VERSION)"' > $@

htslib: 
	$(MAKE) -C htslib

MethIndelRealigner: htslib version.h $(OBJS)
	$(CC) $(OPTS) -I$(INCLUDE_DIRS) $(OBJS) htslib/libhts.a main.c -o MethIndelRealigner -lz -lpthread -lm

clean:
	rm -f *.o MethIndelRealigner

clean-all: clean
	make --directory=htslib clean

install: MethIndelRealigner
	install MethIndelRealigner $(PREFIX)
