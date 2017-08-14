# LATEST MODIFICATION : June 2013
# This is the compiler I will use:
CC = g++
# These are my flags
CFLAGS = -c -O2 -ftree-vectorize #-Wall

oshun.e: state.o formulary.o setup.o vlasov.o export.o main.o 
	$(CC) state.o formulary.o vlasov.o export.o main.o -o oshun.e
	rm main.o

state.o: state.cpp lib-array.h state.h 
	$(CC) $(CFLAGS) state.cpp

formulary.o: formulary.cpp lib-array.h formulary.h 
	$(CC) $(CFLAGS) formulary.cpp

setup.o: setup.cpp lib-array.h lib-algorithms.h state.h formulary.h setup.h
	$(CC) $(CFLAGS) setup.cpp

vlasov.o: vlasov.cpp lib-array.h lib-algorithms.h state.h setup.h vlasov.h 
	$(CC) $(CFLAGS) vlasov.cpp

export.o: export.cpp lib-array.h lib-algorithms.h state.h formulary.h setup.h export.h 
	$(CC) $(CFLAGS) export.cpp

main.o: main.cpp lib-array.h lib-algorithms.h state.h formulary.h setup.h vlasov.h export.h
	$(CC) $(CFLAGS) main.cpp

clean:
	rm -f oshun.e *.o *~ *#
