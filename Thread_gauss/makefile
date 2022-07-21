CFLAGS = -pthread -std=c++14 -O3 -ffast-math
a.out: main.cpp libs/matrix.o libs/thread_function.o
	g++ $(CFLAGS) main.cpp libs/matrix.o libs/thread_function.o -o a.out
libs/matrix.o: libs/matrix.cpp libs/matrix.h
	g++ $(CFLAGS) -c libs/matrix.cpp -o libs/matrix.o
libs/thread_function.o: libs/thread_function.cpp libs/thread_function.h
	g++ $(CFLAGS)  -c libs/thread_function.cpp  -o libs/thread_function.o
clean:
	rm -rf *.o libs/*.o