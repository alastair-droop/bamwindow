CC=gcc
CFLAGS=-Wall
HTSDIR=./htslib

all:
	$(CC) $(CFLAGS) -I$(HTSDIR) -L$(HTSDIR) -o bamwindow bamwindow.c -lhts
	
clean:
	-rm bamwindow
