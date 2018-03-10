# AUTH Raymundo Hernandez Esparza
# DATE March 2017
# REVI March 2017
# REQU no-requirements
#


PRO = cube3D
EXE = $(PRO).x

CC := gcc
CFLAGS := -c -Wunused
LIB := -lm -llapacke -llapack

OBJS := main.o\
	file.o\
        array.o\
	analys.o\
	tableP.o\
	matvec.o\
	lectura.o\
	critical.o\
	mathTools.o\
	pruebaCoef.o\
	numBondPath.o\
	lagrangePol.o

$(PRO) : $(OBJS)
	$(CC) -o $(EXE) $(OBJS) $(LIB)

main.o: main.c file.h array.h analys.h lectura.h critical.h
	$(CC) $(CFLAGS) main.c

file.o: file.c file.h
	$(CC) $(CFLAGS) file.c

array.o: array.c  array.h
	$(CC) $(CFLAGS) array.c

analys.o: analys.c analys.h file.h array.h lectura.h
	$(CC) $(CFLAGS) -DGRAD analys.c

tableP.o: tableP.c tableP.h
	$(CC) $(CFLAGS) tableP.c

matvec.o: matvec.c matvec.h
	$(CC) $(CFLAGS) -DD matvec.c

lectura.o: lectura.c lectura.h file.h array.h 
	$(CC) $(CFLAGS) lectura.c

critical.o: critical.c file.c array.h analys.h critical.h numBondPath.h
	$(CC) $(CFLAGS) critical.c

mathTools.o : mathTools.c mathTools.h
	$(CC) $(CFLAGS) mathTools.c

pruebaCoef.o : pruebaCoef.c pruebaCoef.h
	$(CC) $(CFLAGS) pruebaCoef.c

numBondPath.o: numBondPath.c numBondPath.h file.h array.h
	$(CC) $(CFLAGS) numBondPath.c

lagrangePol.o: lagrangePol.c lagrangePol.h
	$(CC) $(CFLAGS) lagrangePol.c

.PHONY: install
install:
	mv $(EXE) ../bin/.

.PHONY: clean
clean:
	rm $(OBJS) $(EXE)


