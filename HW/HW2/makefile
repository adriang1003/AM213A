FC = gfortran    
FFLAGS = -O3
DFLAGS = -g3
LFLAGS =

OBJECTS = EX13.o
.PHONY: clean

main.ex: $(OBJECTS)
	$(FC) $(LFLAGS) -o main.ex $(OBJECTS)
%.o : %.f90
	$(FC) -c $(DFLAGS) $<
clean:
	rm -f main.ex $(OBJECTS)
	rm -f myplot.dat
	rm -f *.o