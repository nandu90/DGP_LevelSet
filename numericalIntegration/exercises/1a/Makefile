CC = g++
OPT = -o

exercise1a: polylib.o exercise1a.o
	$(CC) $(CPPFLAGS) exercise1a.o polylib.o $(OPT) exercise1a

polylib.o: ../../polylib.c ../../polylib.h
	$(CC) $(CPPFLAGS) -c ../../polylib.c -lm

exercise1a.o: exercise1a.cpp
	$(CC) $(CPPFLAGS) -c exercise1a.cpp -lm


clean:
	rm -rf *.o exercise1a
