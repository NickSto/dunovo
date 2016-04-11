CFLAGS=-Wall

all:
	gcc -Wall -shared -fPIC align.c -o libalign.so
	gcc -Wall -shared -fPIC swalign.c -o libswalign.so -lm
	gcc -Wall -shared -fPIC seqtools.c -o libseqtools.so
	gcc -Wall -shared -fPIC consensus.c -o libconsensus.so

