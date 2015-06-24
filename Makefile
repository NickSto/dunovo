CFLAGS=-Wall

all:
	gcc -Wall -shared -fPIC align.c -o align.so
