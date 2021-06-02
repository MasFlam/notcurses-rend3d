.PHONY: clean

CFLAGS := -O3

rend3d: rend3d.o vecmat.o
	cc -o $@ $^ -lnotcurses-core -lm

rend3d.o: rend3d.c vecmat.o
	cc $(CFLAGS) -c -o $@ rend3d.c

vecmat.o: vecmat.h vecmat.c
	cc $(CFLAGS) -c -o $@ vecmat.c

clean:
	rm -f vecmat.o rend3d.o rend3d
