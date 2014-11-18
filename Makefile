PREFIX = /home/ryand/bin
INCLUDE_DIRS = /home/ryand/include
LIB_DIRS = /home/ryand/lib
CC = gcc
OPTS = -Wall -g #-DDEBUG

OBJS = alignmentHeap.o bloomFilter.o graph.o murmur3.o TargetCreator.o SSW/ssw.o realigner.o SemiGlobal.o

.PHONY: all clean

.SUFFIXES:.c .o

all: TargetCreator Realigner

.c.o:
	$(CC) -c $(OPTS) -I$(INCLUDE_DIRS) $< -o $@

TargetCreator: TargetCreator.o
	$(CC) $(OPTS) -I$(INCLUDE_DIRS) -L$(LIB_DIRS) -o TargetCreator TargetCreator.o TargetCreator_main.c -lhts -lz -lpthread

Realigner: $(OBJS)
	$(CC) $(OPTS) -I$(INCLUDE_DIRS) -L$(LIB_DIRS) $(OBJS) realigner_main.c -o Realigner -lhts -lz -lpthread -lm

clean:
	rm -f *.o TargetCreator Realigner
