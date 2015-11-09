SOURCE= main.cpp polynomial.cpp system.cpp sylvester.cpp vsylvester.cpp
HEADERS= polynomial.h system.h sylvester.h vsylvester.h
OBJS= main.o polynomial.o system.o sylvester.o vsylvester.o
CC= g++
CFLAGS= -c -g
EXEC= equations

all: $(OBJS)
	$(CC) $(OBJS) -o $(EXEC)

main.o: main.cpp system.h sylvester.h vsylvester.h
	$(CC) $(CFLAGS) main.cpp


polynomial.o: polynomial.cpp polynomial.h
	$(CC) $(CFLAGS) polynomial.cpp

system.o: system.cpp system.h polynomial.h
	$(CC) $(CFLAGS) system.cpp

sylvester.o: sylvester.cpp sylvester.h
	$(CC) $(CFLAGS) sylvester.cpp

vsylvester.o: vsylvester.cpp vsylvester.h
	$(CC) $(CFLAGS) vsylvester.cpp


clean:
	rm -rf *o $(EXEC)

count:
	wc $(HEADERS) $(SOURCE)
