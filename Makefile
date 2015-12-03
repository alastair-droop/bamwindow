CC=gcc
CFLAGS=-Wall
HTSDIR=./htslib # Change this to point at your htslib source directory.

bamwindow: htslib
	$(CC) $(CFLAGS) -I$(HTSDIR) -L$(HTSDIR) -o bamwindow bamwindow.c -lhts

htslib:
	cd $(HTSDIR) && $(MAKE)
	
clean:
	-rm bamwindow
	cd $(HTSDIR) && $(MAKE) clean
