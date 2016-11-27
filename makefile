plot: remove
	gnuplot gr.plt
remove: run
	rm prog *.o *.mod
run: comp
	./prog
comp: compmod Prog.f90
	gfortran Prog.f90 DataTypes.o -o prog
compmod: DataTypes.f90
	gfortran -c DataTypes.f90