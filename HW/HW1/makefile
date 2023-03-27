FC = gfortran    
FFLAGS = -O3
DFLAGS = -g3
LFLAGS =

OBJECTS = main.o pi_appx.o
OBJ16 = solbybisection.o bisect.o fcosx.o ftanx.o
.PHONY: clean

all: main.ex output.txt EX12 EX13 EX15 EX16

output.txt: main.ex
	./main.ex > output.txt
main.ex: $(OBJECTS)
	$(FC) $(LFLAGS) -o main.ex $(OBJECTS)
EX12: EX12.o
	$(FC) $(LFLAGS) -o EX12 EX12.o
EX13: EX13.o
	$(FC) $(LFLAGS) -o EX13 EX13.o
EX15: EX15.o
	$(FC) $(LFLAGS) -o EX15 EX15.o
EX16: $(OBJ16)
	$(FC) $(LFLAGS) -o EX16 $(OBJ16)
%.o : %.f90
	$(FC) -c $(DFLAGS) $<
clean:
	rm -f main.ex $(OBJECTS)
	rm -f EX16 $(OBJ16)
	rm -f EX12 EX13 EX15
	rm -f *.o
	rm -f output.txt