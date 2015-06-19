CFLAGS=-Wall

all:
	gcc -shared -fPIC align.c -o align.so
