# AUTH Raymundo Hernandez Esparza
# DATE March 2017
# REVI APRIL 2020
# REQU lapacke / lapack
# VER  v1.2


PRO = cube3D
EXE = $(PRO).x

CC := icc
###CFLAGS := -c -Wunused  -Wno-unused-result -O2 
CFLAGS := -c -Wunused  -O2 
LIB := -lm  -fopenmp

OBJS := main.o\
        file.o\
		  array.o\
		  graph.o\
		  fields.o\
		  transU.o\
		  tableP.o\
		  timing.o\
		  jacobi.o\
		  kernels.o\
		  lagForm.o\
		  lebedev.o\
		  lectura.o\
          findCrit.o\
		  geomData.o\
		  rotation.o\
		  cubeIndex.o\
		  lagrange2.o\
		  replicate.o\
		  mathTools.o\
          numBondPath.o\
          numGradLines.o

$(PRO) : $(OBJS)
	$(CC) $(CFLAGS)  utils.c
	$(CC) -o $(EXE) $(OBJS) utils.o $(LIB)

main.o: main.c file.h array.h  lectura.h struct.h
	$(CC) $(CFLAGS) main.c

file.o: file.c file.h
	$(CC) $(CFLAGS) -DD file.c

array.o: array.c  array.h
	$(CC) $(CFLAGS) array.c

graph.o: graph.c  graph.h
	$(CC) $(CFLAGS) graph.c

fields.o: fields.c fields.h
	$(CC) $(CFLAGS) fields.c -fopenmp

tableP.o: tableP.c tableP.h
	$(CC) $(CFLAGS) tableP.c

timing.o: timing.c timing.h
	$(CC) $(CFLAGS) timing.c

jacobi.o: jacobi.c jacobi.h
	$(CC) $(CFLAGS) jacobi.c

kernels.o: kernels.c kernels.h
	$(CC) $(CFLAGS) -DGRAPH  kernels.c

lagForm.o: lagForm.c lagForm.h
	$(CC) $(CFLAGS) lagForm.c

lebedev.o: lebedev.c lebedev.h
	$(CC) $(CFLAGS) lebedev.c

lectura.o: lectura.c lectura.h
	$(CC) $(CFLAGS) lectura.c

transU.o: transU.c transU.h
	$(CC) $(CFLAGS) transU.c

findCrit.o: findCrit.c findCrit.h
	$(CC) $(CFLAGS)  findCrit.c -fopenmp

critical.o: critical.c critical.h
	$(CC) $(CFLAGS) critical.c

geomData.o: geomData.c geomData.h
	$(CC) $(CFLAGS) geomData.c

rotation.o : rotation.c rotation.h
	$(CC) $(CFLAGS) rotation.c

cubeIndex.o	: cubeIndex.c cubeIndex.h
	$(CC) $(CFLAGS) cubeIndex.c

lagrange2.o : lagrange2.c lagrange2.h
	$(CC) $(CFLAGS) lagrange2.c

replicate.o : replicate.c replicate.h
	$(CC) $(CFLAGS) replicate.c

mathTools.o : mathTools.c mathTools.h
	$(CC) $(CFLAGS) mathTools.c

pruebaCoef.o : pruebaCoef.c pruebaCoef.h
	$(CC) $(CFLAGS) pruebaCoef.c

numBondPath.o: numBondPath.c numBondPath.h 
	$(CC) $(CFLAGS)  numBondPath.c -fopenmp

numGradLines.o: numGradLines.c numGradLines.h
	$(CC) $(CFLAGS)  numGradLines.c -fopenmp

.PHONY: install
install:
	mv $(EXE) $(HOME)/bin/.

.PHONY: clean
clean:
	rm $(OBJS) utils.o $(EXE) *.cube *.dat *.xyz *.log *.pdf *.eps *.gp


