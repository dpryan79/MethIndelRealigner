PREFIX = /home/ryand/bin #This can be changed
INCLUDE_DIRS = htslib
CC = gcc
OPTS = -Wall -g -O3

OBJS = alignmentHeap.o hashtable.o graph.o murmur3.o TargetCreator.o realigner.o SemiGlobal.o threads.o TargetCreator_main.o realigner_main.o
VERSION = 0.1.1HT

#If we're building from a git repo, then append the most recent tag
ifneq "$(wildcard .git)" ""
VERSION := $(VERSION)-$(shell git describe --always --dirty)
endif

.PHONY: all clean htslib install clean-all version.h

.SUFFIXES:.c .o

all: MethIndelRealignerHT

.c.o:
	$(CC) -c $(OPTS) -I$(INCLUDE_DIRS) $< -o $@

version.h:
	echo '#define VERSION "$(VERSION)"' > $@

htslib: 
	$(MAKE) -C htslib

MethIndelRealignerHT: htslib version.h $(OBJS)
	$(CC) $(OPTS) -I$(INCLUDE_DIRS) $(OBJS) htslib/libhts.a main.c -o MethIndelRealignerHT -lz -lpthread -lm

clean:
	rm -f *.o MethIndelRealignerHT

clean-all: clean
	make --directory=htslib clean

install: MethIndelRealignerHT
	install MethIndelRealignerHT $(PREFIX)
