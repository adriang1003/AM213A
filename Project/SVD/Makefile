OBJECTS = Driver_LinAl.o

.PHONY: clean

output.txt: LinAl.exe
	./LinAl.exe > output.txt

LinAl.exe: $(OBJECTS)
	gfortran $(OBJECTS) -llapack -o LinAl.exe

%.o: %.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8  $<

%.mod: %.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8  $<

clean:
	rm -f $(OBJECTS) LinAl.exe
	rm -f output.txt
	rm -f Image_appn_100010.dat
	rm -f Image_appn_100020.dat
	rm -f Image_appn_100040.dat
	rm -f Image_appn_100080.dat
	rm -f Image_appn_100160.dat
	rm -f Image_appn_100320.dat
	rm -f Image_appn_100640.dat
	rm -f Image_appn_101279.dat
	rm -f error.dat