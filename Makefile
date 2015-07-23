CFLAGS=-Wall

all:
	gcc -Wall -shared -fPIC align.c -o align.so
	gcc -Wall -shared -fPIC swalign.c -o swalignc.so -lm
