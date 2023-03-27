OBJECTS = LinAl.o Driver_LinAl.o

MODULES = LinAl.mod

.PHONY: clean

LinAl.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) -o LinAl.exe

%.o: %.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8  $<

%.mod: %.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8  $<

clean:
	rm -f $(OBJECTS) $(MODULES) LinAl.exe
	rm -f linal.mod
	rm -f *.dat