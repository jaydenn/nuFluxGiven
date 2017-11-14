CC = mpiCC
FLAGS = -g -O3 -DMPI -DOMPI_SKIP_MPICXX
LIBS = -lmultinest_mpi -lmpi -llapack -lgfortran -lgsl -lgslcblas -lstdc++ -lm
INCLUDE = -I./src/inc
OBJECTS = src/solFits.o src/fileIO.o src/detectorFunctions.o src/monteCarlo.o src/likelihood.o src/nuRates.o src/globalFit.o src/calcRates.o
default: solFits

solFits: $(OBJECTS)
	$(CC) $(FLAGS) $^ -o $@ $(LIBS)

src/%.o: src/%.cpp
	$(CC) $< $(FLAGS) ${INCLUDE} -c -o $@

clean:
	-rm src/*.o
	-rm -f ./solFits
