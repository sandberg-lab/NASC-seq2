CC = gcc
CFLAGS = -g -Wall

all: hyp1f1_arb_vector.so

hyp1f1_arb_vector.so: hyp1f1_arb_vector.o
	$(CC) hyp1f1_arb_vector.o -shared -o hyp1f1_arb_vector.so -larb -lflint -lmpfr -lgmp

hyp1f1_arb_vector.o: hyp1f1_arb_vector.c
	$(CC) hyp1f1_arb_vector.c -c -fPIC -larb -lflint -lmpfr -lgmp $(CFLAGS) -o hyp1f1_arb_vector.o
