PREFIX = /home/ryand/bin #This can be changed
INCLUDE_DIRS = htslib
CC = gcc
OPTS = -Wall -g #-DDEBUG

OBJS = alignmentHeap.o bloomFilter.o graph.o murmur3.o TargetCreator.o SSW/ssw.o realigner.o SemiGlobal.o

.PHONY: all clean htslib install

.SUFFIXES:.c .o

all: TargetCreator Realigner

.c.o:
	$(CC) -c $(OPTS) -I$(INCLUDE_DIRS) $< -o $@

htslib: 
	$(MAKE) -C htslib

TargetCreator: htslib TargetCreator.o
	$(CC) $(OPTS) -I$(INCLUDE_DIRS) -o TargetCreator TargetCreator.o htslib/libhts.a TargetCreator_main.c -lz -lpthread

Realigner: $(OBJS) htslib
	$(CC) $(OPTS) -I$(INCLUDE_DIRS) $(OBJS) htslib/libhts.a realigner_main.c -o Realigner -lz -lpthread -lm

clean:
	rm -f *.o TargetCreator Realigner
	make --directory=htslib clean

install: TargetCreator Realigner
	install TargetCreator Realigner $(PREFIX)
