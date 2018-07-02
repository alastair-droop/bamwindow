CC=gcc
CFLAGS=-Wall
HTSDIR=/usr/local/lib/

all:
	$(CC) $(CFLAGS) -I$(HTSDIR) -L$(HTSDIR) -o bamwindow bamwindow.c -lm -lhts
	
clean:
	-rm bamwindow
