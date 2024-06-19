CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: es1.exe es2.exe

es1.exe : es1.o random.o
	$(CC) random.o es1.o -o es1.exe

es2.exe : es2.o random.o
	$(CC) random.o es2.o -o es2.exe
	

es1.o : es1.cpp
	$(CC) -c es1.cpp -o es1.o $(CFLAGS)

es2.o : es2.cpp
	$(CC) -c es2.cpp -o es2.o $(CFLAGS)
	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o es1.exe es2.exe seed.out