module task
  use mpi
  implicit none
contains
  subroutine GetMaxCoordinates(A, x1, y1, x2, y2)
    implicit none
    real(8), intent(in), dimension(:,:) :: A
    integer(4), intent(out) :: x1, y1, x2, y2
    integer(4) :: n, L, R, Up, Down, m, tmp
    real(8), allocatable :: current_column(:), B(:,:)
    real(8) :: current_sum, max_sum,max_sum_globvalue
    logical :: transpos
    integer(4) :: mpiErr, mpiSize, mpiRank , tmp_mpiRank, global_mpiRank! MPI
    
    call mpi_comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)
    
    call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)
    
    m = size(A, dim=1) 
    n = size(A, dim=2) 
    transpos = .FALSE.
    
    if (m < n) then 
       transpos = .TRUE.   
       B = transpose(A)
       m = size(B, 1) 
       n = size(B, 2) 
    else
       B = A     
    endif
    
    allocate(current_column(m))
    
    
    max_sum=B(1,1)
    x1=1
    y1=1
    x2=1
    y2=1
    
    !        do L=1, n
    do L = mpiRank + 1, n, mpiSize ! Параллелим цикл do в соответствии с количеством процессов  
       current_column = B(:, L)            
       do R=L,n
          
          if (R > L) then 
             current_column = current_column + B(:, R)
          endif
          
          call FindMaxInArray(current_column, current_sum, Up, Down) 
          
          if (current_sum > max_sum) then
             max_sum = current_sum
             x1 = Up
             x2 = Down
             y1 = L
             y2 = R
          endif
                   
       end do
    end do
        
    !   Делаем свёртку по максимальному элементу
    call mpi_reduce(max_sum, max_sum_globvalue, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, mpiErr)
    
    !   Рассылаем максимальный элемент всем процессам
    call mpi_bcast(max_sum_globvalue, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpiErr)
    
    !   Проверяем, из какого процесса мы взяли максимальную сумму
    tmp_mpiRank = -1
    if (max_sum .eq. max_sum_globvalue) then
       tmp_mpiRank = mpiRank
    endif
    call mpi_allreduce(tmp_mpiRank, global_mpiRank, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpiErr)
    
    call mpi_bcast(x1,1, MPI_REAL8, global_mpiRank, MPI_COMM_WORLD, mpiErr)  
    call mpi_bcast(x2,1, MPI_REAL8, global_mpiRank, MPI_COMM_WORLD, mpiErr)
    call mpi_bcast(y1,1, MPI_REAL8, global_mpiRank, MPI_COMM_WORLD, mpiErr)
    call mpi_bcast(y2,1, MPI_REAL8, global_mpiRank, MPI_COMM_WORLD, mpiErr)
    
    deallocate(current_column)  
    
    if (transpos) then  
       tmp = x1
       x1 = y1
       y1 = tmp
       
       tmp = y2
       y2 = x2
       x2 = tmp
    endif
    
  end subroutine GetMaxCoordinates
  
  
  subroutine FindMaxInArray(a, Sum, Up, Down)
    real(8), intent(in), dimension(:) :: a
    integer(4), intent(out) :: Up, Down
    real(8), intent(out) :: Sum
    real(8) :: cur_sum
    integer(4) :: minus_pos, i
    !integer(4) :: mpiErr, mpiSize, mpiRank ! MPI
    
    Sum = a(1)
    Up = 1
    Down = 1
    cur_sum = 0
    minus_pos = 0
    
    
    
    do i=1, size(a)
       cur_sum = cur_sum + a(i)
       if (cur_sum > Sum) then
          Sum = cur_sum
          Up = minus_pos + 1
          Down = i
       endif
       
       if (cur_sum < 0) then
          cur_sum = 0
          minus_pos = i
       endif
       
    enddo

  end subroutine FindMaxInArray

  
end module task



