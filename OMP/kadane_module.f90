module kadane_module
  use omp_lib
  implicit none
contains
  subroutine GetMaxCoordinates(A, x1, y1, x2, y2)
    implicit none
    real(8), intent(in), dimension(:,:) :: A
    integer(4), intent(out) :: x1, y1, x2, y2
    integer(4) :: n, L, R, Up, Down, m, tmp
    real(8), allocatable :: current_column(:), B(:,:)
    real(8) :: current_sum, max_sum
    logical :: transpos
    
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
    
!    call omp_set_num_threads(4)

     !$omp parallel shared(x1, y1, x2, y2, max_sum, B, n) default(private)
     !$omp do schedule(dynamic)
    
    do L=1, n
       
       current_column = B(:, L)            
       do R=L,n
          
          if (R > L) then 
             current_column = current_column + B(:, R)
          endif
          
          call FindMaxInArray(current_column, current_sum, Up, Down) 
          
          if (current_sum > max_sum) then
          
          !$omp critical         
          if (current_sum > max_sum) then
             max_sum = current_sum
             x1 = Up
             x2 = Down
             y1 = L
             y2 = R
          endif
          !$omp end critical
          endif
          
       end do
    end do
    
    !$omp end do
    !$omp end parallel
    
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
    integer(4) :: mpiErr, mpiSize, mpiRank ! MPI

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

  
end module kadane_module



