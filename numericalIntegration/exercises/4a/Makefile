CC = g++
OPT = -o

code: polylib.o code.o
	$(CC) $(CPPFLAGS) code.o polylib.o $(OPT) code

polylib.o: ../../polylib.c ../../polylib.h
	$(CC) $(CPPFLAGS) -c ../../polylib.c -lm

code.o: code.cpp
	$(CC) $(CPPFLAGS) -c code.cpp -lm


clean:
	rm -rf *.o code
