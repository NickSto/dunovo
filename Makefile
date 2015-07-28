CFLAGS=-Wall

all:
	gcc -Wall -shared -fPIC alignc.c -o alignc.so
	gcc -Wall -shared -fPIC swalignc.c -o swalignc.so -lm
	gcc -Wall -shared -fPIC seqtoolsc.c -o seqtoolsc.so
	gcc -Wall -shared -fPIC consensusc.c -o consensusc.so

