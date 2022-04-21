objects = predP.o combP.o omegaMap.o io.o mcmc.o tools.o analyse.o likelihood.o
myutils = argumentwizard.h controlwizard.h DNA.h lotri_matrix.h \
	matrix.h MLST.h mutation.h mydouble.h myerror.h myutils.h \
	random.h utils.h vector.h lissim.h
headers = omegaMap.h $(myutils)
idir = "./"
OLEVEL = 0

## INSTRUCT MAKE TO MAKE MULTIPLE GOALS
.PHONY : all
all : predictiveP combineP

## LINK
predictiveP : $(objects)
	gcc -w -O3 -o predictiveP $(objects) -lstdc++

combineP : $(objects)
	gcc -w -O3 -o combineP $(objects) -lstdc++

## COMPILE
combP.o : combineP.cpp $(headers)
	gcc -w -O3 -c -o combP.o combineP.cpp -I$(idir)

predP.o : predP.cpp $(headers)
	gcc -w -O3 -c -o predP.o predP.cpp -I$(idir)

omegaMap.o : omegaMap.cpp $(headers)
	gcc -w -O3 -c -o omegaMap.o omegaMap.cpp -I$(idir)

io.o : io.cpp $(headers)
	gcc -w -O3 -c -o io.o io.cpp -I$(idir)

likelihood.o : likelihood.cpp $(headers)
	gcc -w -O3 -c -o likelihood.o likelihood.cpp -I$(idir)

mcmc.o : mcmc.cpp $(headers)
	gcc -w -O$(OLEVEL) -c -o mcmc.o mcmc.cpp -I$(idir)

tools.o : tools.c paml.h
	gcc -w -O3 -c -o tools.o tools.c -I$(idir)

## MAKE CLEAN
clean :
	-rm $(objects)
