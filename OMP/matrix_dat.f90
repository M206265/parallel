module matrix_dat
contains
  subroutine make_data(n)
    integer i, j, n
    real tmp
    real, allocatable, dimension(:,:) :: matrix
    open(10, file = 'data1.dat')
!    write(10,"('# ',I6)") n
    write(10, *) n
    allocate(matrix(n,n))
    call RANDOM_NUMBER(matrix)
    matrix = matrix*10 - 5
    write(10, *) matrix
    close(10)
    deallocate(matrix)

  end subroutine make_data
end module matrix_dat
