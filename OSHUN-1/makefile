# LATEST MODIFICATION : July 28th 2011
# This is the compiler I will use:
CC = mpic++
# These are my flags
CFLAGS = -c -O3

OSHUN-2D-im.e: def-Input.o def-Nmethods.o def-State.o def-Setup.o def-Output.o def-Export.o def-Parallel.o def-VlasovMax.o def-Actions.o def-RK.o def-ImplicitE.o def-FokkerPlanck.o def-PLSource.o main.o 
	$(CC) def-Input.o def-Nmethods.o def-State.o def-Setup.o def-Output.o def-Export.o def-Parallel.o def-VlasovMax.o def-Actions.o def-RK.o def-ImplicitE.o def-FokkerPlanck.o def-PLSource.o main.o -o OSHUN-2D-im.e
	rm def-Input.o
	rm def-Setup.o
	rm def-Output.o
	rm main.o

def-Input.o: def-Input.cpp matrices.h decl-input.h 
	$(CC) $(CFLAGS) def-Input.cpp

def-Nmethods.o: def-Nmethods.cpp matrices.h decl-nmethods.h 
	$(CC) $(CFLAGS) def-Nmethods.cpp

def-State.o: def-State.cpp matrices.h decl-input.h decl-state.h 
	$(CC) $(CFLAGS) def-State.cpp

def-Setup.o: def-Setup.cpp matrices.h decl-input.h decl-state.h decl-setup.h 
	$(CC) $(CFLAGS) def-Setup.cpp

def-Output.o: def-Output.cpp matrices.h decl-input.h decl-state.h decl-output.h 
	$(CC) $(CFLAGS) def-Output.cpp

def-Export.o: def-Export.cpp matrices.h decl-input.h decl-state.h decl-export.h 
	$(CC) $(CFLAGS) def-Export.cpp

def-Parallel.o: def-Parallel.cpp matrices.h decl-input.h decl-state.h decl-output.h decl-export.h decl-parallel.h
	$(CC) $(CFLAGS) def-Parallel.cpp

def-VlasovMax.o: def-VlasovMax.cpp matrices.h decl-input.h decl-state.h decl-vlasovmax.h
	$(CC) $(CFLAGS) def-VlasovMax.cpp

def-PLSource.o: def-PLSource.cpp matrices.h decl-input.h decl-state.h decl-plsource.h
	$(CC) $(CFLAGS) def-PLSource.cpp

def-Actions.o: def-Actions.cpp matrices.h decl-input.h decl-state.h decl-vlasovmax.h decl-actions.h
	$(CC) $(CFLAGS) def-Actions.cpp

def-RK.o: def-RK.cpp matrices.h decl-input.h decl-state.h decl-vlasovmax.h decl-actions.h decl-RK.h
	$(CC) $(CFLAGS) def-RK.cpp

def-FokkerPlanck.o: def-FokkerPlanck.cpp matrices.h decl-nmethods.h decl-input.h decl-state.h decl-fokkerplanck.h
	$(CC) $(CFLAGS) def-FokkerPlanck.cpp

def-ImplicitE.o: def-ImplicitE.cpp matrices.h decl-input.h decl-state.h decl-vlasovmax.h decl-actions.h decl-RK.h decl-fokkerplanck.h decl-implicitE.h 
	$(CC) $(CFLAGS) def-ImplicitE.cpp

main.o: main.cpp matrices.h decl-input.h decl-state.h decl-setup.h decl-export.h decl-parallel.h decl-vlasovmax.h decl-plsource.h decl-actions.h decl-RK.h decl-implicitE.h decl-fokkerplanck.h 
	$(CC) $(CFLAGS) main.cpp

clean:
	rm -f OSHUN-2D-im.e *.o *~ *#

