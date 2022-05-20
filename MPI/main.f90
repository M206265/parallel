program main
  use task
  use matrix_dat
  implicit none
  real(8), allocatable, dimension(:,:) :: matrix_A
  integer(4) :: n, x1, x2, y1, y2
  real(8) :: start, finish 
  integer(4) :: mpiErr


  
!  call make_data(2000)
  
  open(1, file='data1.dat')
  read(1, *) n
  allocate(matrix_A(n,n))
  read(1,*) matrix_A
  close(1)
call mpi_init(mpiErr) ! mpi
  call GetMaxCoordinates(matrix_A, x1, y1, x2, y2)
call mpi_finalize(mpiErr) ! MPI
  write(*,*) x1, y1, x2, y2

  deallocate(matrix_A)

  
end program main

