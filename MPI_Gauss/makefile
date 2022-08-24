CFLAGS = -O3
a.out: libs/MPI_matrix.o libs/matrix.o main.cpp
	mpic++ $(CFLAGS) libs/MPI_matrix.o libs/matrix.o main.cpp -o a.out
libs/MPI_matrix.o: libs/MPI_matrix.cpp libs/matrix.o libs/MPI_matrix.h
	mpic++ $(CFLAGS) libs/matrix.o -c libs/MPI_matrix.cpp -o libs/MPI_matrix.o
libs/matrix.o: libs/matrix.cpp libs/matrix.h
	g++ $(CFLAGS) -c libs/matrix.cpp -o libs/matrix.o
clean:
	rm -rf *.o libs/*.o