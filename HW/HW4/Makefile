OBJECTS = LinAl.o Driver_LinAl.o

MODULES = LinAl.mod

.PHONY: clean

output.txt: LinAl.exe
	./LinAl.exe > output.txt

LinAl.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) -o LinAl.exe

%.o: %.f90
	gfortran -c -freal-8-real-4  $<

%.mod: %.f90
	gfortran -c -freal-8-real-4  $<

clean:
	rm -f $(OBJECTS) $(MODULES) LinAl.exe
	rm -f output.txt
	rm -f linal.mod
	rm -f ChoDegree3.dat
	rm -f ChoDegree5.dat
	rm -f QRDegree3.dat
	rm -f QRDegree5.dat
