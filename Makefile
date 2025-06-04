all:
	mpic++ $(F).cpp -o $(F)

sec:
	g++ $(F).cpp -o $(F)

run4:
	mpirun -np 4 ./$(F)

run8:
	mpirun -np 8 ./$(F)

clean:
	rm -f ./$(F)

omp:
	g++ -fopenmp $(F).cpp -o $(F)

ompi:
	mpic++ -fopenmp $(F).cpp -o $(F)