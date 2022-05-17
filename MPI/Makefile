files=$(wildcard *.f90)
obj=$(patsubst %.f90, %.o, $(files))

res: result clean

result: main
	time ./main > result 
main: $(obj)                                                          
	mpif90 $^ -o $@
%.o: %.f90
	mpif90 -c $^ -o $@ 
clean:
	rm -f *.o *.mod


