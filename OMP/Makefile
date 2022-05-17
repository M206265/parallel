files=$(wildcard *.f90)
obj=$(patsubst %.f90, %.o, $(files))

res: result clean

result: main
	time ./main > result 
main: $(obj)                                                          
	gfortran $^ -o $@ -fopenmp
%.o: %.f90
	gfortran -c $^ -o $@
#main.o: kadane_module.mod matrix_dat.mod
#kadane_module.mod: omp_lib.mod
clean:
	rm -f *.o *.mod


